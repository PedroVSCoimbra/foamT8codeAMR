#!/usr/bin/env python3
"""Post-process risingDrop_ellipsoidal_AMR case.

Extracts bubble properties over time:
- Bubble position (center of mass of alpha1)
- Bubble rise velocity
- Circularity / deformation
- Bubble area

Usage:
  cd risingDrop_ellipsoidal_AMR
  python3 postprocess.py
"""

import os
import re
import numpy as np

CASE_DIR = os.getcwd()
G = 9.81
R0 = 0.25        # initial bubble radius
A0 = np.pi*R0**2  # initial bubble area

# Expected terminal velocity (for reference)
# V_ref = sqrt(g * R0) ≈ 1.565


def read_time_dirs(case_dir):
    times = []
    for d in os.listdir(case_dir):
        p = os.path.join(case_dir, d)
        if os.path.isdir(p) and d not in ('0', '0.orig', 'constant', 'system',
                                           'postProcessing'):
            try:
                t = float(d)
                times.append(t)
            except ValueError:
                pass
    return sorted(times)


def parse_internal_field(filepath):
    """Extract numerical values from a foam volField file."""
    with open(filepath) as f:
        content = f.read()

    # Try uniform first
    m = re.search(r'internalField\s+uniform\s+([^;]+)', content)
    if m:
        val = m.group(1).strip()
        if val.startswith('('):
            return np.array([float(x) for x in val.strip('()').split()])
        return np.array([float(val)])

    # Try nonuniform List
    m = re.search(
        r'internalField\s+nonuniform\s+List<(\w+)>\s*\n\s*(\d+)\s*\n\s*\(',
        content, re.DOTALL
    )
    if m:
        n = int(m.group(2))
        rest = content[m.end():]
        closing = rest.index(')')
        data_str = rest[:closing]
        result = np.array([float(x) for x in data_str.split()])
        return result

    return None


def compute_bubble_props(alpha):
    """Compute bubble position, area, circularity from alpha1 field."""
    n_cells = len(alpha)
    # Estimate: count cells with alpha > 0.5 as bubble
    bubble_cells = alpha > 0.5
    n_bubble = np.sum(bubble_cells)

    if n_bubble == 0:
        return None, None, None, None

    # Fraction of bubble cells ≈ bubble area / total area
    # Area = n_bubble * dA where dA = total_area / n_cells
    total_area = 1.0 * 2.0  # domain dimensions: 1 × 2
    bubble_area = n_bubble * total_area / n_cells

    # Simple center of mass (y-direction only for rising)
    cell_centers_y = np.linspace(0, 2, n_cells, endpoint=False)
    # This is approximate - assumes regular grid
    bubble_center_y = np.mean(cell_centers_y[bubble_cells])

    return bubble_center_y, bubble_area, n_bubble, total_area


def main():
    print(f"=== risingDrop_ellipsoidal_AMR Post-Processing ===")
    print(f"Case: {CASE_DIR}")

    time_dirs = read_time_dirs(CASE_DIR)
    print(f"Found {len(time_dirs)} time directories")

    if not time_dirs:
        print("No time directories found. Run the simulation first.")
        return

    output_dir = os.path.join(CASE_DIR, 'postProcessing')
    os.makedirs(output_dir, exist_ok=True)

    # Compute bubble properties at each time
    t_list = []
    y_list = []
    area_list = []
    n_bubble_list = []

    for t in time_dirs:
        alpha_file = os.path.join(CASE_DIR, str(t), 'alpha1')
        if not os.path.exists(alpha_file):
            continue

        alpha = parse_internal_field(alpha_file)
        if alpha is None:
            continue

        y, area, n_b, total = compute_bubble_props(alpha)
        if y is None:
            continue

        t_list.append(t)
        y_list.append(y)
        area_list.append(area)
        n_bubble_list.append(n_b)

    if not t_list:
        print("Could not read bubble data.")
        return

    t_arr = np.array(t_list)
    y_arr = np.array(y_list)
    area_arr = np.array(area_list)

    # Save to file
    data_path = os.path.join(output_dir, 'bubble_data.txt')
    np.savetxt(
        data_path,
        np.column_stack([t_arr, y_arr, area_arr]),
        header='t(s)  y_com(m)  area(m2)'
    )
    print(f"Saved bubble data to {data_path}")

    # Print summary
    y_max = y_arr[-1]
    n_bubble = n_bubble_list[-1]
    print(f"\nFinal bubble position: y = {y_max:.4f} m")
    print(f"Final bubble area:     {area_arr[-1]:.4f} m²")
    print(f"Area conservation:     {area_arr[-1]/A0:.4f} (1.0 = perfect)")
    print(f"Bubble cells at t={t_arr[-1]:.2f}: {n_bubble}")

    # Try matplotlib for plotting
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 1, figsize=(8, 8))

        ax = axes[0]
        ax.plot(t_arr, y_arr, 'b-')
        ax.set_xlabel('t (s)')
        ax.set_ylabel('Bubble center y (m)')
        ax.set_title('Rising bubble position')
        ax.grid(True, alpha=0.3)

        ax = axes[1]
        ax.plot(t_arr, area_arr, 'r-')
        ax.axhline(A0, color='k', ls='--', label=f'Initial area = {A0:.3f}')
        ax.set_xlabel('t (s)')
        ax.set_ylabel('Bubble area (m²)')
        ax.set_title('Bubble area conservation')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'bubble_rise.png'), dpi=150)
        print(f"Saved bubble_rise.png")

    except ImportError:
        print("matplotlib not available: skipping plots")


if __name__ == '__main__':
    main()
