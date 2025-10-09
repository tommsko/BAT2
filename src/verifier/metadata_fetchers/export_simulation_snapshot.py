import os
import sys
import subprocess
import tempfile
import shutil


def generate_snapshot(tpr_path, gro_path, output_path):
    if not os.path.exists(tpr_path):
        raise FileNotFoundError(f"TPR file not found: {tpr_path}")

    if gro_path and not os.path.exists(gro_path):
        raise FileNotFoundError(f"GRO file not found: {gro_path}")

    with tempfile.TemporaryDirectory() as tmpdir:
        vmd_script = os.path.join(tmpdir, "render.tcl")
        image_path = os.path.join(tmpdir, "snapshot.tga")

        # VMD Tcl script
        with open(vmd_script, "w") as f:
            f.write(f"""
mol new "{gro_path or tpr_path}" type {{gro}} waitfor all
mol addfile "{tpr_path}" type {{tpr}} waitfor all

# Select and split water and solute
set water [atomselect top "resname SOL or water"]
set not_water [atomselect top "not resname SOL and not water"]

$not_water set beta 0.0
$water set beta 1.0

$not_water drawframes all all
$water drawframes all all

# Representations
mol delrep 0 top
mol representation VDW 1.0 12.0
mol color Name
mol selection "not resname SOL and not water"
mol material Transparent
mol addrep top

mol representation VDW 0.6 12.0
mol color Name
mol selection "resname SOL or water"
mol material Ghost
mol addrep top

# Center and zoom manually
set sel [atomselect top "all"]
set center [measure center $sel weight mass]
set minmax [measure minmax $sel]
set vecmin [lindex $minmax 0]
set vecmax [lindex $minmax 1]

# Calculate bounding box size
set size [vecsub $vecmax $vecmin]
set maxdim [lindex [lsort -real -decreasing $size] 0]

# Move molecule to center
$sel moveby [vecscale -1 $center]

# Reset view and apply rotation of 45 degrees around Y axis
display resetview

# Rotate 45 degrees around Y axis
set angle_rad [expr {{75 * acos(-1) / 180.0}}]
set c [expr {{cos($angle_rad)}}]
set s [expr {{sin($angle_rad)}}]

# Build rotation matrix around Y axis
set rot_mat {{
    $c 0 $s 0
    0 1 0 0
    -$s 0 $c 0
    0 0 0 1
}}

display matrix $rot_mat

# Zoom out: scale so molecule fits smaller in bigger frame (2x bigger)
set zoomscale [expr {{0.5 / $maxdim}}]
scale $zoomscale

# Display settings
display resize 800 800
display projection orthographic
display depthcue off
display antialias on
color Display Background white

# Render
render TachyonInternal {image_path}

exit
""")


        # Run VMD
        try:
            subprocess.run(["vmd", "-dispdev", "text", "-e", vmd_script], check=True)
        except subprocess.CalledProcessError as e:
            print("VMD failed:", e)
            sys.exit(1)

        # Convert TGA to PNG (to ensure transparency)
        from PIL import Image
        with Image.open(image_path) as img:
            img = img.convert("RGBA")
            img.save(output_path)

