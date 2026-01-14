import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# --------------------------------------------------
# Page setup
# --------------------------------------------------
st.set_page_config(page_title="GC-MS vs Integrated System Simulation", layout="wide")
st.title("Simulation: Integrated Analytical System vs Traditional GC-MS")

# --------------------------------------------------
# User notice
# --------------------------------------------------
st.warning(
    "⚠️ هذه المحاكاة تعليمية: تعطي تصورًا عن الأجهزة وتحليلها، "
    "وليست تجربة حقيقية في المختبر."
)

# --------------------------------------------------
# System selection
# --------------------------------------------------
system_choice = st.radio(
    "Select system to simulate:",
    ("Integrated System (Your Device)", "Traditional GC-MS")
)

# --------------------------------------------------
# Sample input
# --------------------------------------------------
st.header("Sample Input")

compound_names = st.text_area(
    "Enter compound names (comma separated):",
    "Acetone, Methanol, Ethanol"
)

boiling_points = st.text_area(
    "Enter boiling points °C (same order):",
    "56, 65, 78"
)

ion_eff = st.text_area(
    "Enter ionization efficiency % (same order):",
    "70, 60, 65"
)

# Clean inputs
compound_names = [c.strip() for c in compound_names.split(",") if c.strip()]
boiling_points = [float(b.strip()) for b in boiling_points.split(",") if b.strip()]
ion_eff = [float(i.strip())/100 for i in ion_eff.split(",") if i.strip()]

if not (len(compound_names) == len(boiling_points) == len(ion_eff)):
    st.error("⚠️ Please enter the same number of values for all fields.")
    st.stop()

n_compounds = len(compound_names)

# --------------------------------------------------
# Thermal & system settings
# --------------------------------------------------
if system_choice == "Integrated System (Your Device)":
    st.header("Sample Preparation & Conditioning (Integrated System)")

    # Sample prep oven (does NOT affect chromatogram directly)
    prep_temp = st.slider("Sample Prep Oven Temperature (°C)", 50, 400, 150)
    prep_pressure = st.slider("Sample Prep Oven Pressure (kPa)", 10.0, 500.0, 101.3)

    # Hot channel & µTD
    hot_channel_temp = st.slider("Hot Channel / µTD Temperature (°C)", 100, 400, 250)

    column_rate = st.slider("Column heating rate (°C/min)", 1, 50, 5)

    # Hot-factor: improves focusing & transfer efficiency
    hot_factor = 1 + (hot_channel_temp - 150)/400 + (column_rate/50)*0.1

else:
    st.header("Traditional GC-MS Settings")
    column_rate = st.slider("Column heating rate (°C/min)", 1, 50, 5)

# --------------------------------------------------
# Time axis
# --------------------------------------------------
x = np.linspace(0, 20, 1200)

# --------------------------------------------------
# Signal containers
# --------------------------------------------------
y_column_ion = np.zeros_like(x)
y_column_neutral = np.zeros_like(x)
peak_table = []

# --------------------------------------------------
# Resolution function
# --------------------------------------------------
def resolution(t1, w1, t2, w2):
    return 2 * abs(t2 - t1) / (w1 + w2)

# --------------------------------------------------
# Simulation loop
# --------------------------------------------------
for i in range(n_compounds):
    # Retention time (approx.)
    # Depends mainly on boiling point, could include polarity, column length, pressure
    center = 2 + i*2 + (boiling_points[i]/100)

    # Peak width (approx.)
    # Linked to boiling point and column heating rate
    width = 0.15 + (boiling_points[i]/500) + (10/column_rate)*0.02

    if system_choice == "Integrated System (Your Device)":
        # Narrower & stronger peaks due to Hot Channel / µTD
        width_column = width / hot_factor
        intensity_ion = ion_eff[i] * hot_factor
        intensity_neutral = (1 - ion_eff[i]) * hot_factor
    else:
        width_column = width
        intensity_ion = ion_eff[i]
        intensity_neutral = 0

    # Gaussian peaks
    peak_ion = intensity_ion * np.exp(-((x - center)**2)/(2*width_column**2))
    peak_neutral = intensity_neutral * np.exp(-((x - center)**2)/(2*width_column**2))

    y_column_ion += peak_ion
    y_column_neutral += peak_neutral

    peak_table.append({
        "Compound": compound_names[i],
        "Retention Time (min)": round(center,2),
        "Peak Width": round(width_column,3),
        "Ionic Intensity": round(intensity_ion,3),
        "Neutral Intensity": round(intensity_neutral,3)
    })

# --------------------------------------------------
# Resolution calculation
# --------------------------------------------------
for i in range(1, len(peak_table)):
    t1 = peak_table[i-1]["Retention Time (min)"]
    t2 = peak_table[i]["Retention Time (min)"]
    w1 = peak_table[i-1]["Peak Width"]
    w2 = peak_table[i]["Peak Width"]
    peak_table[i]["Resolution vs Previous"] = round(resolution(t1,w1,t2,w2),2)

peak_table[0]["Resolution vs Previous"] = None

df = pd.DataFrame(peak_table)

# --------------------------------------------------
# Plot
# --------------------------------------------------
st.header("Simulated Chromatogram")

fig, ax = plt.subplots(figsize=(14,6))

# Show column signal
ax.plot(x, y_column_ion, color="blue", label="Column – Ionic Path")
if system_choice == "Integrated System (Your Device)":
    ax.plot(x, y_column_neutral, color="green", linestyle="--", label="Column – Neutral Path (Integrated Device)")

# Highlight µTD / Hot Channel importance
if system_choice == "Integrated System (Your Device)":
    for row in peak_table:
        ax.axvspan(row["Retention Time (min)"]-0.1, row["Retention Time (min)"]+0.1,
                   color="orange", alpha=0.1, label="µTD / Hot Channel Effect" if row==peak_table[0] else "")

# Annotate compounds
for row in peak_table:
    ax.scatter(row["Retention Time (min)"], row["Ionic Intensity"], zorder=3, color="blue")
    ax.text(row["Retention Time (min)"], row["Ionic Intensity"]+0.02, row["Compound"], ha="center")

ax.set_xlabel("Retention Time (min)")
ax.set_ylabel("Signal (a.u.)")
ax.set_title(system_choice)
ax.legend()
st.pyplot(fig)

# --------------------------------------------------
# Results table
# --------------------------------------------------
st.header("Detected Peaks & Resolution")
st.dataframe(df)

# --------------------------------------------------
# Download CSV
# --------------------------------------------------
csv = df.to_csv(index=False).encode('utf-8')
st.download_button(
    label="Download peaks table as CSV",
    data=csv,
    file_name='simulated_peaks.csv',
    mime='text/csv'
)

# --------------------------------------------------
# Scientific notes
# --------------------------------------------------
st.info(
    "🔹 This simulation is educational. Hot Channel and µTD are crucial for sample focusing "
    "and safe transfer to the column, improving separation and sample integrity, "
    "but do NOT directly add extra peaks. Retention time, peak width, and intensity are approximated; "
    "real GC-MS behavior depends on column chemistry, length, carrier gas, pressure, and detector physics."
)
