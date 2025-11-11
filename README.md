## Supplementary Information for `DM_Gomp.ini` (DarkHistory/CLASS Integration)

This branch, `DH_CLASS_integration`, includes Dark Matter decay integration using the DarkHistory.

To successfully run the included parameter file, `DM_Gomp.ini`, **without a full installation of the DarkHistory code**, please follow these steps to simulate the expected data environment.

-----

### 1\. Set up the Data Directory Structure

Ensure the `DarkHistory` directory is **at the same level** as the `class_gomp` directory. You will need to create the required subdirectories within it:

```bash
# Navigate to the directory containing both class_gomp and DarkHistory
# If DarkHistory doesn't exist yet, this creates the necessary structure:
mkdir -p DarkHistory/data
mkdir -p DarkHistory/decay_e
```

-----

### 2\. Place the Data File

Move the required pre-calculated DH data file into the `DarkHistory/decay_e` folder.

```bash
mv class_gomp/1.000000_24.000000_67.660000_0.260690_0.048970_0.060000_6.000000_-2.040000.txt DarkHistory/decay_e/
```

-----

### 3\. Set the Environment Variable

You must set the `DH_DATA_DIR` environment variable to point to the base data directory you just created.

```bash
# NOTE: Update the path with the correct absolute path on the collaborator's machine.
export DH_DATA_DIR=/path/to/file/DarkHistory/data
```

-----

### 4\. Run CLASS

You can now run the `DM_Gomp.ini` parameter file from the `class_gomp` directory:

```bash
# Navigate to your class_gomp directory
cd class_gomp/

# Then run the code
./class DM_Gomp.ini
```
