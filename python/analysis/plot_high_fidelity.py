#!/usr/bin/env python3
"""M24 high-fidelity environment plots: ground track + gravity/ephemeris."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 9,
    "figure.titlesize": 12,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
})


def main():
    root = Path(__file__).resolve().parents[2]
    data = root / "data"
    figs = root / "artifacts" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(data / "m24_environment.csv")

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    fig.suptitle("M24 — Ground Track + Gravity Model Comparison (100 min LEO)", fontweight="bold")
    axs[0].plot(df["lon_deg"], df["lat_deg"], color="royalblue")
    axs[0].set_title("Geodetic Ground Track (ECEF)")
    axs[0].set_xlabel("Longitude (deg)")
    axs[0].set_ylabel("Latitude (deg)")
    axs[0].grid(True)
    axs[1].plot(df["time_s"] / 60.0, df["j234_total_mag"] - 8.425, color="forestgreen",
                label="J2/J3/J4 total minus 8.425 baseline")
    axs[1].set_title("Total Gravity (detail, central removed for display)")
    axs[1].set_xlabel("Time (min)")
    axs[1].set_ylabel("m/s^2 above 8.425")
    axs[1].grid(True)
    axs[1].legend()
    fig.tight_layout()
    fig.savefig(figs / "m24_ground_track.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
    fig.suptitle("M24 — Fixed vs Ephemeris Moon (100 min)", fontweight="bold")
    ax.plot(df["time_s"] / 60.0, df["moon_fixed_err_km"], color="crimson")
    ax.set_title("Fixed-Moon Error Growth")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Error (km)")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(figs / "m24_ephemeris.png")
    plt.close(fig)
    print("Wrote M24 figures")


if __name__ == "__main__":
    main()
