import os
import statistics
import tomli
from concurrent.futures import ProcessPoolExecutor, as_completed
import math
from tqdm import tqdm
 

# =========================================================
# CONFIG
# =========================================================

with open("config.toml", "rb") as f:
    CFG = tomli.load(f)

BASE = CFG["paths"]["base_dir"]
ODB_DIR = os.path.join(BASE, CFG["paths"]["odb_dir"])
ST_DIR = os.path.join(BASE, CFG["paths"]["st_dir"])
META_FILE = os.path.join(BASE, CFG["paths"]["metadata"])

STATS_FILE = CFG["files"]["stats"]
WLIST_FILE = CFG["files"]["whitelist"]

RUN_ODB = CFG["processing"]["run_odb_extraction"]
MAX_PAR = CFG["processing"]["max_parallel"]
PREF_LEN = CFG["processing"]["station_prefix_len"]
BAD_STD = CFG["qc"]["invalid_std"]


# =========================================================
# HELPERS
# =========================================================

def read_single_line_file(filename, skip_first=False):
    rows = []
    with open(filename) as f:
        if skip_first:
            next(f, None)
        for line in f:
            if line.strip():
                rows.append(line.split())
    return rows

def station4(s):
    return s.strip().strip("'\"")[:4]

def station8(s):
    return s.strip().strip("'\"")[:8]

def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0  # Earth radius in km

    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)

    a = (math.sin(dphi / 2) ** 2 +
         math.cos(phi1) * math.cos(phi2) *
         math.sin(dlambda / 2) ** 2)

    return 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a))


# =========================================================
# STAGE 0 (OPTIONAL): ODB ? *.st
# =========================================================

def extract_station_rows(file_path, station_id):
    rows = read_single_line_file(file_path, skip_first=True)
    out = []
    for r in rows:
        if r[0] == station_id:
            out.append((r[0], float(r[4]), float(r[5])))
    return out


def write_st_file(station, rows):
    with open(os.path.join(ST_DIR, f"{station}.st"), "a") as f:
        for sid, v1, v2 in rows:
            f.write(f"{sid:<8} {v1:9.4f} {v2:10.6f}\n")

def process_single_odb(args):
    file_path, station_id = args
    rows = extract_station_rows(file_path, station_id)
    if rows:
        write_st_file(station_id[1:-1], rows)


def run_odb_extraction():
    print("Running ODB -> *.st extraction")

    stations = read_single_line_file(META_FILE)
    station_ids = [s[0] for s in stations]

    if CFG["files"]["odb_extract_prefix"]:
       ODB_PREFIX=CFG["files"]["odb_extract_prefix"]
    else:
       ODB_PREFIX="ecma_"

    odb_files = [
        os.path.join(ODB_DIR, f)
        for f in os.listdir(ODB_DIR)
        if f.lower().startswith(ODB_PREFIX)
    ]

    tasks = [
        (odb, sid)
        for sid in station_ids
        for odb in odb_files
    ]

    total = len(tasks)

    with ProcessPoolExecutor(max_workers=MAX_PAR) as exe:
        futures = [exe.submit(process_single_odb, t) for t in tasks]

        for _ in tqdm(
            as_completed(futures),
            total=total,
            desc="ODB extraction",
            unit="task"
        ):
            pass


        for f in as_completed(futures):
            f.result()


# =========================================================
# STAGE 1: *.st ? stats
# =========================================================

def compute_stats():
    """
    Reads all *.st files in ST_DIR and computes
    mean / std / skew per processing center.

    The station key is derived from the CONTENT of the file
    (first column), normalized to 8 characters:
        <station(4)><network(4)>
    """

    stats = {}

    # reset stats file
    with open(STATS_FILE, "w"):
        pass

    for fname in sorted(os.listdir(ST_DIR)):
        if not fname.lower().endswith(".st"):
            continue

        path = os.path.join(ST_DIR, fname)
        rows = read_single_line_file(path)

        if not rows:
            continue

        # Station ID from data 
        station_id8 = station8(rows[0][0])
        station_id4 = station4(rows[0][0])
        
         

        values = []
        for r in rows:
            try:
                values.append(float(r[2]))
            except (IndexError, ValueError):
                continue

        if not values:
            continue

        mean = statistics.mean(values) * 1000
        median = statistics.median(values) * 1000
        std = statistics.stdev(values) * 1000 if len(values) > 1 else BAD_STD
        skew = 3 * (mean - median) / std if std < BAD_STD else 0.0

        stats[fname] = {
            "station4": station_id4,
            "station8": station_id8,
            "mean": mean,
            "std": std,
            "skew": skew,
        }

        with open(STATS_FILE, "a") as f:
            f.write(
                f"{fname:<8} {mean:9.4f} {std:10.6f} {skew:10.6f}\n"
            )

    return stats

# =========================================================
# STAGE 2: stats + metadata ? wlist
# =========================================================

def normalize_station_id(s):
    return s.strip().strip("'\"")

def load_station_metadata():
    meta = {}

    for row in read_single_line_file(META_FILE):
        if len(row) < 4:
            continue

        sid = normalize_station_id(row[0])

        meta[sid] = {
            "lat": float(row[1]),
            "lon": float(row[2]),
            "alt": float(row[3]),
        }

    return meta

def apply_min_distance_filter(candidates, meta, min_km):
    """
    candidates: dict[station4] -> stats dict
    meta: station8 -> {lat, lon, alt}

    Returns filtered candidates.
    """

    accepted = []

    # Sort by quality FIRST (best stations get priority)
    ordered = sorted(
        candidates.values(),
        key=lambda s: s["std"]
    )

    for s in ordered:
        st8 = s["station8"]
        m = meta.get(st8)
        if not m:
            continue

        lat, lon = m["lat"], m["lon"]

        too_close = False
        for acc in accepted:
            m2 = meta[acc["station8"]]
            d = haversine_km(lat, lon, m2["lat"], m2["lon"])
            if d < min_km:
                too_close = True
                break

        if not too_close:
            accepted.append(s)

    return {s["station4"]: s for s in accepted}


def write_whitelist(stats, meta):
    """
    One whitelist entry per PHYSICAL station (4 chars),
    selecting the best processing center (8 chars).
    """

    with open(WLIST_FILE, "w"):
        pass

    best = {}

    # 1. select best processing center per physical station
    for fname, s in stats.items():
        if s["std"] >= BAD_STD:
            continue

        st4 = s["station4"]

        if st4 not in best or s["std"] < best[st4]["std"]:
            best[st4] = s

    if CFG["spatial_filter"]["enable"]:
    	best = apply_min_distance_filter(
        	best,
        	meta,
        	CFG["spatial_filter"]["min_distance_km"]
   	)

    # 2. write whitelist
    for st4, s in sorted(best.items()):
        st8 = s["station8"]

        if st8 not in meta:
            print(f"WARNING: missing metadata for {st8}")
            continue

        m = meta[st8]

        with open(WLIST_FILE, "a") as f:
            f.write(
                f"{st8:<8} "
                f"{m['lat']:6.2f} "
                f"{m['lon']:6.2f} "
                f"{m['alt']:7.2f} "
                f"{15.0:4.1f} "
                f"{s['mean']/1000:8.4f} "
                f"{s['std']:7.2f} "
                f"{0.1:4.2f}\n"
            )

# =========================================================
# MAIN
# =========================================================

if __name__ == "__main__":

    os.chdir(BASE)
    os.makedirs(ST_DIR, exist_ok=True)

    if RUN_ODB:
        run_odb_extraction()
    else:
        print("Skipping ODB extraction (using existing *.st files)")

    stats = compute_stats()
 
    ## optionally, outprint the numebr of duplciates for each station
    # from collections import Counter
    # print("Duplicate stations in stats:",
    #  {k: v for k, v in Counter(s["station4"] for s in stats.values()).items() if v > 1})

    meta = load_station_metadata()
    write_whitelist(stats, meta)

    print("Whitelist computation finished successfully.")

