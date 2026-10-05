#!/usr/bin/env python3
import os, csv, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 14,
    "mathtext.fontset": "stix",
    "text.usetex": True,
})

CASE_DIR = Path.cwd()
G = 9.81; H0 = 0.6; a = 0.6

def _load_csv(path):
    t, h = [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            t.append(float(row['t_s'])); h.append(float(row['h_m']))
    return np.array(t), np.array(h)

exp_time_h1, exp_value_h1 = _load_csv(str(CASE_DIR/'data_x1_2725m.csv'))
exp_time_h2, exp_value_h2 = _load_csv(str(CASE_DIR/'data_x2_2228m.csv'))
t_max = max(exp_time_h1[-1], exp_time_h2[-1])

probe_dir = CASE_DIR/'postProcessing'/'probes'/'0'
probe_p = None
if probe_dir.exists() and (probe_dir/'p').exists():
    raw = np.loadtxt(probe_dir/'p', skiprows=1)
    if raw.ndim > 1: probe_p = (raw[:,0], raw[:,2])

def parse_fvmatrix(content):
    m = re.search(r'internalField\s+nonuniform\s+List<\w+>\s*\n\s*(\d+)\s*\n\s*\(([\s\S]*?)\)\s*;', content)
    if not m:
        um = re.search(r'internalField\s+uniform\s+([^;]+)', content)
        if um:
            v = um.group(1).strip()
            if v.startswith('('):
                vs = [float(x) for x in v.strip('()').split()]
                return np.array([np.linalg.norm(vs)]) if len(vs) > 1 else np.array([vs[0]])
            return np.array([float(v)])
        return None
    vals = []
    for line in m.group(2).strip().split('\n'):
        line = line.strip()
        if not line: continue
        if line.startswith('('):
            vv = [float(x) for x in line.strip('()').split()]
            vals.append(np.linalg.norm(vv) if len(vv) > 1 else vv[0])
        else:
            for v in line.split():
                try: vals.append(float(v))
                except ValueError: pass
    return np.array(vals) if vals else None

def parse_ascii_list(fp):
    with open(fp) as f: c = f.read()
    c = re.sub(r'/\*.*?\*/', '', c, flags=re.DOTALL)
    c = re.sub(r'//.*', '', c)
    idx = c.rfind('}')
    if idx >= 0: c = c[idx+1:]
    c = c.strip()
    m = re.match(r'(\d+)\s*\(', c)
    if not m: return None
    n_exp = int(m.group(1))
    start = m.end() - 1
    end = c.rfind(')')
    if end <= start: return None
    body = c[start+1:end].strip()

    face_items = re.findall(r'\d+\([^)]*\)', body)
    if len(face_items) == n_exp:
        return [[int(x) for x in item.split('(')[1].strip(')').split()] for item in face_items]

    vec_items = re.findall(r'\([^)]*\)', body)
    if len(vec_items) == n_exp:
        return [[float(x) for x in item.strip('()').split()] for item in vec_items]

    vals = re.findall(r'[-+]?\d+', body)
    if len(vals) == n_exp: return [[int(v)] for v in vals]

    vals2 = body.split()
    if len(vals2) == n_exp:
        try: return [[int(v)] for v in vals2]
        except: pass
        try: return [[float(v)] for v in vals2]
        except: pass
    return None

def read_cell_centers(td):
    poly = CASE_DIR/td/'polyMesh'
    if not poly.exists(): poly = CASE_DIR/'constant'/'polyMesh'
    if not poly.exists(): return None
    pts = parse_ascii_list(poly/'points')
    faces = parse_ascii_list(poly/'faces')
    owner = parse_ascii_list(poly/'owner')
    if pts is None or faces is None or owner is None: return None
    pts_np = np.array([p[:3] for p in pts], dtype=float)
    oi = [o[0] for o in owner]
    n_cells = max(oi)+1
    neigh = parse_ascii_list(poly/'neighbour')
    if neigh is not None: n_cells = max(n_cells, max(n[0] for n in neigh)+1)
    cs = np.zeros((n_cells,3)); cw = np.zeros(n_cells, dtype=int)
    for fi, fp in enumerate(faces):
        ci = oi[fi]; fc = np.mean(pts_np[fp], axis=0)
        cs[ci] += fc; cw[ci] += 1
    if neigh is not None:
        for fi, ni in enumerate(neigh):
            ci = ni[0]; fc = np.mean(pts_np[faces[fi]], axis=0)
            cs[ci] += fc; cw[ci] += 1
    valid = cw > 0
    cc = np.full((n_cells,3), np.nan)
    cc[valid] = cs[valid] / cw[valid, np.newaxis]
    return cc

def water_heights_amr(td, probe_xs=(2.725, 2.228)):
    af = CASE_DIR/td/'alpha1'
    if not af.exists(): return None, None
    with open(af) as f: a = parse_fvmatrix(f.read())
    if a is None: return None, None
    cc = read_cell_centers(td)
    if cc is None: return None, None
    valid = ~np.isnan(cc[:,0])
    cc = cc[valid]
    a = a[valid]
    results = []
    tol = 0.05
    for px in probe_xs:
        strip = np.where(np.abs(cc[:,0]-px) < tol)[0]
        if len(strip) == 0: results.append(0.0)
        else:
            above = a[strip] > 0.5
            if not np.any(above): results.append(0.0)
            else: results.append(float(np.max(cc[strip[above], 1])))
    return tuple(results)

def read_time_dirs():
    return sorted([d.name for d in CASE_DIR.iterdir()
        if d.is_dir() and d.name not in ('0','0.orig','constant','system','postProcessing','__pycache__')
        and not d.name.startswith('.')], key=float)

def plot_water_levels(times, h1v, h2v, out):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
    ax1.plot(exp_time_h1, exp_value_h1, 'o', label=r'Zhainakov \& Kurbanaliev (2013)',
             markersize=8, markerfacecolor='none', markeredgecolor='black',
             markeredgewidth=1.5, alpha=0.7)
    if times and h1v: ax1.plot(times, h1v, '-', label='OpenFOAM')
    ax1.set_xlabel('$t$ [s]'); ax1.set_ylabel('$h$ [m]')
    ax1.legend(fontsize=10); ax1.grid(True, alpha=0.3)

    ax2.plot(exp_time_h2, exp_value_h2, 'o', label=r'Zhainakov \& Kurbanaliev (2013)',
             markersize=8, markerfacecolor='none', markeredgecolor='black',
             markeredgewidth=1.5, alpha=0.7)
    if times and h2v: ax2.plot(times, h2v, '-', label='OpenFOAM')
    ax2.set_xlabel('$t$ [s]'); ax2.set_ylabel('$h$ [m]')
    ax2.legend(fontsize=10); ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out/'water_levels.png', dpi=150)
    print('Saved water_levels.png')

def main():
    print(f'=== damBreakValidation AMR ===\nCase: {CASE_DIR}')
    out = CASE_DIR/'postProcessing'; out.mkdir(exist_ok=True)

    tds = read_time_dirs()
    times, h1v, h2v = [], [], []
    for td in tds:
        t = float(td)
        if t > t_max: break
        h1, h2 = water_heights_amr(td)
        if h1 is not None:
            times.append(t); h1v.append(h1); h2v.append(h2)
        if len(times) % 20 == 0 and len(times) > 0:
            print(f'  processed {len(times)}/{len(tds)} time dirs...')

    print(f'Processed {len(times)} time directories (up to t={t_max:.2f}s)')
    if probe_p is not None:
        print(f'Probe pressure data: {len(probe_p[0])} points')
    plot_water_levels(times, h1v, h2v, out)

if __name__ == '__main__':
    main()
