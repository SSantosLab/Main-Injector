#this is a script for determining OneRing's performance
#based on the math in Sec. 2.4 of the ongoing MI paper 
#written by Elise Fuller in 2026 w/ some help from Claude Opus 5.5 (Claude mostly wrote docstrings and comments, like the one below, and helped with debugging. the docstrings, comments, and debugging were manually checked and changed where needed)
"""
onering_metrics.py; goal is to test performance of OneRing-generated observing plan
 
Implements the math of the "Offline testing and validation" section of the
MI paper:
 
    pi_h   = integral of the skymap prob density over hex tile T_h (hex.prob in the hex class)
    
    a_i    = |F_i| / |T_h|          fractional coverage of dither i. F_i is the region of T_h covered by DECam's footprint at dither i. 

    
    a_12   = |F_1 union F_2| / |T_h|    region imaged on BOTH passes
 
    P_dis  = sum_{h in H} pi_h * a_12 * min(d1, d2)             H is set of hexes observed on both passes. this is the prob counterpart is detected in at least two epochs for obs plan with two dithers per hex (2 passes)
    P_dis  = sum_h pi_h * [ min(d1,d2) a_h^(1,1) + d2 a_h^(2+) ]          (>2 passes) this works if we implement a plan with >2 passes 
    P_plan = P_dis with every hex observed on all passes, no time limit (ideal scenario)
 
A point of the tile is "discovered" if it is imaged in >= 2 epochs, and under
the paper's nested-depth assumption ("a source detectable in the shallow epoch
is always detectable in the deeper one") the probability of that is the
2nd-largest d among the epochs that cover that point. 
 
"""

from __future__ import annotations
 
import numpy as np
import pandas as pd
from matplotlib.path import Path
from scipy.spatial import cKDTree

#reproduced list of DECam dithers
dithers_list = [
    [0.0, 0.0],
    [0.06389, 0.287436],
    [0.0484, -0.6725],
    [-0.2395875, -0.135264],
    [0.9423775, -0.405792],
    [-0.76668, 0.473424],
    [-0.543065, -0.828492],
]


# Paper table, cumulative coverage after k dithers (for calibration checks).
PAPER_TOTAL_COVERAGE = [0.8522, 0.9636, 0.9857, 0.9932, 0.9977, 0.9990, 0.9995]


ACTIVE_CCDS = [c for c in range(1, 63) if c not in (2, 61)]
 
 
def dither_index(dither, tol=1e-4):
    """Map a dither offset [d_ra, d_dec] back to its index in dithers_list (-1 if unknown)."""
    d = np.asarray(dither, float)
    for k, ref in enumerate(dithers_list):
        if np.allclose(d, ref, atol=tol):
            return k
    return -1


def read_ccd_corners(path="zeroed_corners.list"):
    """Read zeroed CCD corners -> {ccd: (4,2) array of [off_ra, off_dec] in deg}."""
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    out = {}
    for _, row in df.iterrows():
        ccd = int(float(row["CCD"]))
        out[ccd] = np.array([[row[f"Corner {k} Ra"], row[f"Corner {k} Dec"]] for k in (1, 2, 3, 4)], float)
    return out

#unit vecs for nearest-neighbor calculations on sphere
def _unit_xyz(ra, dec):
    ra, dec = np.deg2rad(np.asarray(ra, float)), np.deg2rad(np.asarray(dec, float))
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


class FootprintGeometry:
    """
    Pixel masks of each dithered footprint F_i restricted to the hex tile T_h.
 
    Parameters
    ----------
    ccd_polygons : list of (4,2) arrays, footprint of the working CCDs at dither [0,0]
                   (offsets in degrees, flat-sky).
    tile_polygon : (N,2) array defining T_h. (that's a hex in an observing plan!) 
    dithers      : list of [d_ra, d_dec] offsets 
    n_pix        : raster resolution along the long axis of the tile. this is set to be the same as OneRing hexalation
 
    Use `calibration_report()` to compare a_1 and cumulative coverage with the
    paper table... if a_1 is not ~0.852 your T_h differs from the one
    used for the paper's dither table.
    """
 
    def __init__(self, ccd_polygons, tile_polygon, dithers=dithers_list, n_pix=700):
        self.dithers = [np.asarray(d, float) for d in dithers]

        self.tile_polygon = np.asarray(tile_polygon, float)
 
        xmin, ymin = self.tile_polygon.min(0)
        xmax, ymax = self.tile_polygon.max(0)
        step = max(xmax - xmin, ymax - ymin) / n_pix
        xs = np.arange(xmin + step / 2, xmax, step)
        ys = np.arange(ymin + step / 2, ymax, step)
        X, Y = np.meshgrid(xs, ys)
        pts = np.column_stack([X.ravel(), Y.ravel()])
        in_tile = Path(self.tile_polygon).contains_points(pts)
        self.pts = pts[in_tile]              # pixels of T_h (equal area)
        self.n_tile_pix = len(self.pts)
 
        paths = [Path(np.asarray(p, float)) for p in ccd_polygons]
        self.masks = np.zeros((len(self.dithers), self.n_tile_pix), bool)
        for k, d in enumerate(self.dithers):
            shifted = self.pts - d           # point in tile covered if (pt - dither) in footprint
            m = np.zeros(self.n_tile_pix, bool)
            for p in paths:
                m |= p.contains_points(shifted)
            self.masks[k] = m
 
    @classmethod
    def from_corners_file(cls, path="zeroed_corners.list", tile="outline", **kw):
        #this gets the outline of a hex w/out dithers
        corners = read_ccd_corners(path)
        active = [corners[c] for c in ACTIVE_CCDS if c in corners]
        if tile == "outline":# T_h = camera outline decam2hp sums over
            import decam2hp
            ora, odec = decam2hp.cameraOutlineAtZero("decam")
            return cls(active, tile_polygon=np.column_stack([ora, odec]), **kw)
        return cls(active, tile_polygon=np.asarray(tile, float), **kw)  # explicit polygon

    def restricted_to(self, keep):
        """used for per-hex T_h."""
        g = object.__new__(FootprintGeometry)
        g.dithers, g.tile_polygon = self.dithers, self.tile_polygon
        g.pts, g.masks = self.pts[keep], self.masks[:, keep]
        g.n_tile_pix = int(keep.sum())
        return g
        
    #a values from paper
    def a(self, i):
        """a_i = |F_i| / |T_h|"""
        return self.masks[i].mean()
 
    def a_overlap(self, i, j):
        """a_ij = |F_i union F_j| / |T_h|   (a_12 = a_overlap(0, 1))"""
        return (self.masks[i] & self.masks[j]).mean()
 
    def a_union(self, idx):
        return np.any(self.masks[list(idx)], axis=0).mean()
 
    def a11_a2plus(self, n_passes):
        """(a^(1,1), a^(2+)) for passes = dithers 0..n_passes-1, first pass = dither 0."""
        n1 = self.masks[0]
        n2 = self.masks[1:n_passes].sum(0) if n_passes > 1 else np.zeros_like(n1, int)
        return (n1 & (n2 == 1)).mean(), (n2 >= 2).mean()
 
    def calibration_report(self):
        rows = []
        for k in range(len(self.dithers)):
            rows.append(dict(dither=k, offset=list(self.dithers[k]), a_k=self.a(k),
                             cumulative=self.a_union(range(k + 1)),
                             paper_cumulative=PAPER_TOTAL_COVERAGE[k] if k < len(PAPER_TOTAL_COVERAGE) else np.nan))
        return pd.DataFrame(rows)
 
    def discovery_fraction(self, epochs, n_required=2):
        """
        Fraction of T_h (weighted by detection prob.) imaged in >= n_required epochs.
 
        epochs : list of (dither_index, d) for every exposure of this hex.
        Returns mean over tile pixels of the n_required-th largest d among the
        epochs covering that pixel (0 if fewer than n_required cover it).
        n_required=2 is the paper's P_dis definition; n_required=1 gives the
        single-epoch (d-weighted) coverage.
        """
        if len(epochs) < n_required:
            return 0.0
        idx = np.array([e[0] for e in epochs], int)
        d = np.array([e[1] for e in epochs], float)
        if np.any(idx < 0):
            raise ValueError(f"unknown dither in epochs {epochs}")
        vals = np.where(self.masks[idx], d[:, None], -np.inf)    # (E, P)
        kth = -np.sort(-vals, axis=0)[n_required - 1]            # n-th largest per pixel
        kth[~np.isfinite(kth)] = 0.0
        return float(kth.mean())
 
    def covered_fraction(self, dither_indices):
        """Plain localisation coverage (union) of a tile."""
        if len(dither_indices) == 0:
            return 0.0
        return self.a_union(dither_indices)
 
 
def _epochs_by_hex(plan, until_mjd=None, strict_end=True):
    """plan DataFrame -> {hex_id: [(dither_index, d), ...]} in time order."""
    p = plan
    if until_mjd is not None:
        col = "mjd_end" if strict_end else "mjd_start"
        p = p[p[col] <= until_mjd]
    p = p.sort_values("mjd_start")
    out = {}
    for hid, di, d in zip(p["hex_id"], p["dither_index"], p["detP"]):
        out.setdefault(hid, []).append((int(di), float(d)))
    return out

def decam2hp_geometries(hexes, corners_path, hexFile, n_pix=400, search_deg=3.0):
    #old way of  finding geometry w/ same method as OneRing
    #slightly different, allows hex_id: geometry to differ per hex
    
    import decam2hp 
    import hexalate
    corners = read_ccd_corners(corners_path)
    active = [corners[c] for c in ACTIVE_CCDS if c in corners]
    raC, decC, _ = hexalate.getHexCenters(hexFile)
    xyzC = _unit_xyz(raC, decC)
    treeC = cKDTree(xyzC)
    ora, odec = decam2hp.cameraOutlineAtZero("decam")
    # CCD masks are identical for every hex in the offset frame -> compute once on a region
    # slightly larger than the camera outline, then cut out each hex's own T_h
    base = FootprintGeometry(active, tile_polygon=np.column_stack([ora, odec]) * 1.05, n_pix=n_pix)
    chord = 2 * np.sin(np.deg2rad(search_deg) / 2)

    geoms = {}
    for hid, ra_h, dec_h in zip(hexes["hex_id"], hexes["ra"], hexes["dec"]):
        k = treeC.query(_unit_xyz([ra_h], [dec_h]))[1][0]        
        rk, dk = raC[k], decC[k]
        # offsets in the CCD-corner frame -> sky, same convention as place_ccd()
        sky = np.column_stack([rk + base.pts[:, 0] / np.cos(np.deg2rad(dk)), dk + base.pts[:, 1]])
        keep = decam2hp.hexPath(rk, dk, "decam").contains_points(sky)
        for j in treeC.query_ball_point(xyzC[k], chord):
            if j < k:                                              
                keep &= ~decam2hp.hexPath(raC[j], decC[j], "decam").contains_points(sky)
        geoms[hid] = base.restricted_to(keep)
    return geoms


def _geom_for(geom, hid):
    """geom is either one FootprintGeometry for all hexes or a {hex_id: geometry} dict."""
    return geom[hid] if isinstance(geom, dict) else geom
def p_dis(plan, hexes, geom, until_mjd=None, n_required=2, strict_end=True):
    """
    P_dis of an executed plan (paper eq.).  `plan` = OneRing plan log,
    `hexes` = OneRing hex table (needs hex_id, prob).  Exposures ending after
    until_mjd are dropped (strict_end=True) -- use it to emulate a shorter night.
    """
    pi = hexes.set_index("hex_id")["prob"]
    ep = _epochs_by_hex(plan, until_mjd, strict_end)
    return float(sum(pi[h] * _geom_for(geom, h).discovery_fraction(e, n_required) for h, e in ep.items()))
 
 
def p_loc(plan, hexes, geom, until_mjd=None, strict_end=True):
    """Localisation probability imaged at least once (no detection prob.)."""
    pi = hexes.set_index("hex_id")["prob"]
    ep = _epochs_by_hex(plan, until_mjd, strict_end)
    return float(sum(pi[h] * _geom_for(geom, h).covered_fraction([e[0] for e in es]) for h, es in ep.items()))
 
 
def p_plan(hexes, geom, n_passes=2):
    """P_plan: every hex in `hexes` observed with dithers 0..n_passes-1, no time limit."""
    total = 0.0
    for hid, prob, d1, d2 in zip(hexes["hex_id"], hexes["prob"], hexes["d1"], hexes["d2"]):
        epochs = [(0, d1)] + [(k, d2) for k in range(1, n_passes)]
        total += prob * _geom_for(geom, hid).discovery_fraction(epochs)
    return float(total)
 
 
# Closed forms of p_dis for two passes and 2+ passes. likely will only use two passes, but abilitiy to do more is maintained here
def p_dis_two_pass_closed_form(hexes_observed_twice, geom):
    a12 = geom.a_overlap(0, 1)
    return float(sum(r.prob * a12 * min(r.d1, r.d2) for r in hexes_observed_twice.itertuples()))
 
 
def p_dis_multi_pass_closed_form(hexes, geom, n_passes):
    a11, a2p = geom.a11_a2plus(n_passes)
    return float(sum(r.prob * (min(r.d1, r.d2) * a11 + r.d2 * a2p) for r in hexes.itertuples()))
 
 
def p_dis_curve(plan, hexes, geom, t0=None):
    """Cumulative P_dis, P_loc after each exposure ends (step function in time)."""
    cols = ["t_hours", "mjd", "P_dis", "P_loc", "n_exp"]
    if len(plan) == 0:
        return pd.DataFrame(columns=cols)
    pi = hexes.set_index("hex_id")["prob"]
    p = plan.sort_values("mjd_end")
    t0 = p["mjd_start"].min() if t0 is None else t0
    epochs, c_dis, c_loc = {}, {}, {}
    tot_dis = tot_loc = 0.0
    rows = [(0.0, t0, 0.0, 0.0, 0)]
    for n, (hid, di, d, te) in enumerate(zip(p["hex_id"], p["dither_index"], p["detP"], p["mjd_end"]), 1):
        epochs.setdefault(hid, []).append((int(di), float(d)))
        g = _geom_for(geom, hid)
        new_dis = pi[hid] * g.discovery_fraction(epochs[hid])        # only this hex changed
        new_loc = pi[hid] * g.covered_fraction([e[0] for e in epochs[hid]])
        tot_dis += new_dis - c_dis.get(hid, 0.0)
        tot_loc += new_loc - c_loc.get(hid, 0.0)
        c_dis[hid], c_loc[hid] = new_dis, new_loc
        rows.append(((te - t0) * 24, te, tot_dis, tot_loc, n))
    return pd.DataFrame(rows, columns=cols)
 

def validate_plan(plan, hexes, min_secondpass_wait_s=180.0, max_airmass=1.8,
                  stoptime=None, tol_s=1.0):
    """Return a dict of rule violations / diagnostics. Empty lists == all good."""
    issues = {}
    p = plan.sort_values("mjd_start").reset_index(drop=True)
    gap = (p["mjd_start"].values[1:] - p["mjd_end"].values[:-1]) * 86400
    issues["overlapping_exposures"] = list(np.where(gap < -tol_s)[0])
 
    h = hexes.set_index("hex_id")
    vis_bad, wait_bad, dither_bad = [], [], []
    for hid, g in p.groupby("hex_id"):
        g = g.sort_values("mjd_start")
        if "rise_mjd" in h and "set_mjd" in h:
            bad = (g["mjd_start"] < h.at[hid, "rise_mjd"] - tol_s / 86400) | \
                  (g["mjd_start"] > h.at[hid, "set_mjd"] + tol_s / 86400)
            vis_bad += list(g.index[bad])
        if len(g) >= 2:
            w = (g["mjd_start"].values[1] - g["mjd_end"].values[0]) * 86400
            if w < min_secondpass_wait_s - tol_s:
                wait_bad.append((hid, w))
        if list(g["dither_index"]) != list(range(len(g))):
            dither_bad.append((hid, list(g["dither_index"])))
    issues["exposure_outside_visibility"] = vis_bad
    issues["second_pass_too_soon"] = wait_bad
    issues["dither_order_wrong"] = dither_bad
    if "airmass" in p and p["airmass"].notna().any():
        issues["airmass_too_high"] = list(p.index[p["airmass"] > max_airmass])
    if stoptime is not None:
        issues["exposures_ending_after_stop"] = list(p.index[p["mjd_end"] > stoptime + tol_s / 86400])
 
    n_epochs = p.groupby("hex_id").size()
    orphans = n_epochs[n_epochs == 1].index
    issues["n_orphan_first_passes"] = int(len(orphans))
    issues["prob_in_orphan_first_passes"] = float(h.loc[orphans, "prob"].sum()) if len(orphans) else 0.0
    return issues
 
 
def summarise_run(plan, hexes, geom, n_passes=2, stoptime=None, t0=None):
    #probably the best thing to use for testing!!
    #plan is a pandas dataframe containing an observing plan, as the code is written rn. 
    #.          required columns: "hex_id, dither_index, detP, exptime, mjd_start, mjd_end, slew time
    #within onering, you'd incorporate this with, e.g,
    #plan_rows.append(dict(hex_id=best_hex.hex_id, dither_index=M.dither_index(best_hex.dither),
                              # detP=best_hex.detP[p], exptime=exptime,
                              # mjd_start=current_mjd, mjd_end=new_mjd, slew_time=slewtime))
    #appending each row as a hex is added to the observing plan.

    #other inputs: 
    #hexes: a dataframe of all possible hexes on the skymap; it has columns for hex_id, prob, d1 (first pass detection prob from strategy), d2 (2nd pass det prob from strategy), rise_mjd (hex rise time), set_mjd (hex set time), and region (inner or outer). sample line defining it:
    # hexes_df = pd.DataFrame([dict(hex_id=h.hex_id, prob=h.prob, d1=h.detP[0], d2=h.detP[1],
    #                           rise_mjd=h.rise_mjd, set_mjd=h.set_mjd,
    #                           region='inner' if h in og_inner_hexlist else 'outer')
    #                      for h in all_hexes])
    
    #geom: this defines the DECam footprint itself. you can compute this quickly for OneRing with the FootprintGeometry.from_corners_file function here with 'zeroed_corners.list'. 
    #stoptime: time of sunrise
    #t0: start time for an observing night


    #count how many exposures each hex got (count # of rows with each hexid)
    #result is series indexed vy hexid
    n_epochs = plan.groupby("hex_id").size() if len(plan) else pd.Series(dtype=int)

    #P_dis of the plan as scheduled. strictend=True means end time is manually set
    pd_ = p_dis(plan, hexes, geom, until_mjd=stoptime, strict_end=stoptime is not None)

    #P_plan; what we'd get if every hex got every pass and no time limit (ideal conditions!)
    pp = p_plan(hexes, geom, n_passes)

    #separate inner (high prob) and outer (low prob) hexes 
    inner = hexes[hexes["region"] == "inner"] if "region" in hexes else hexes.iloc[:0]

    #use first exposure start time if didn't manually add start time
    t0 = plan["mjd_start"].min() if (t0 is None and len(plan)) else t0
    return dict(
        #num of hexes labelled inner vs outer in the plan
        n_hexes_inner=int((hexes.get("region") == "inner").sum()) if "region" in hexes else np.nan,
        n_hexes_outer=int((hexes.get("region") == "outer").sum()) if "region" in hexes else np.nan,
        #total number of planned exposures
        n_exposures=int(len(plan)),
        #total number of hexes observed at least once
        n_hexes_observed=int(len(n_epochs)),
        #total number of hexes observed 2+ times. these are the only ones that count toward P_dis
        n_hexes_2plus=int((n_epochs >= 2).sum()),
        #check what fraction of inner hexes were observed at least twice
        frac_inner_2plus=float(np.mean([n_epochs.get(h, 0) >= 2 for h in inner["hex_id"]])) if len(inner) else np.nan,
        #total skymap probability in every hex created by hexalation, whether it was observed or not
        #in other words... the prob density integrated over the whole skymap
        sum_pi_plan=float(hexes["prob"].sum()),
        #total skymap probability of just hexes observed 2+ times
        sum_pi_observed_2plus=float(hexes.set_index("hex_id").loc[n_epochs[n_epochs >= 2].index, "prob"].sum()) if len(n_epochs) else 0.0,
        #prob the counterpart is imaged at least once; ignores detection probability
        P_loc=p_loc(plan, hexes, geom, until_mjd=stoptime, strict_end=stoptime is not None),
        #calculated above. P_dis of plan itself
        P_dis=pd_,
        #this is P_dis if we didn't require the counterpart is imaged at least twice
        P_dis_1ep=p_dis(plan, hexes, geom, until_mjd=stoptime, n_required=1, strict_end=stoptime is not None),
        #found above; ideal conditions
        P_plan=pp,
        #THIS IS THE BEST METRIC!!!! evaluates how well plan went compared to ideal conditions/unlimited time 
        efficiency=pd_ / pp if pp > 0 else np.nan,

        #total amount of observing exposure time (open shutter time)
        open_shutter_h=float(plan["exptime"].sum() / 3600) if len(plan) else 0.0,
        #total time the plan took in reality 
        wall_clock_h=float((plan["mjd_end"].max() - t0) * 24) if len(plan) else 0.0,
        #how much time did we spend slewing?
        total_slew_s=float(plan["slew_time"].sum()) if "slew_time" in plan and len(plan) else np.nan,

        #these r the constants from the paper and defined at the top of this code
        #fraction of hex covered by one pointing w/out dither
        a1=np.mean([g.a(0) for g in geom.values()]) if isinstance(geom, dict) else geom.a(0),
        #frction of hex covered by two passes
        a12=np.mean([g.a_overlap(0, 1) for g in geom.values()]) if isinstance(geom, dict) else geom.a_overlap(0, 1),
    )
