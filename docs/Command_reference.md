# Reference for extract_isosurface_height_map_parallel.py

```
positional arguments:
  infile_afm            AFM frequency shift file at different heigts, produced with ppafm-plot_results --save_df
  infile_height         Height profile file for a specific iso surface, generated with extract_isosurface_height_map.py
  outfile               Name for the constant iso AFM image outfile

optional arguments:
  -h, --help            show this help message and exit
  --offset OFFSET, -o OFFSET
                        constant offset in angstrom at which the constant iso AFM image will be calculated. If the constant height AFM images are adjusted to fit the charge density, this is not necessary. Can be useful to find the correct parameters for the AFM calculation.
  --no_stm              supresses the output of the corresponding height profile (constant current STM image)
  --print_effective_height, -peh
                        saves height map used in the AFM image (can differ from image of the --no_stm flag due to the possibly different cells of AFM and STM calculation and due to interpolation)
  --print_effective_index, -pei
                        same as --print_effective_height, but outputs the used grid indices instead of the actual heights. Used mainly for errorshooting.
  --no_freq_interpolation, -nfi
                        disable frequency interpolation between heights. Used mainly for errorshooting.
  --no_height_interpolation, -nhi
                        disable height interpolation from height profile to AFM grid. Used mainly for errorshooting.
  --full_folder, -ff    searches for all height files a folder (includes subdirectories). Folder is given by infile_height argument. Overwrites outfile
  --afm_as_txt, -aat    prints the AFM image as a txt file with the format row column frequency_shift
  --stm_as_txt, -sat    prints the STM image as a txt file with the format row column height. Works only if --no_stm is not set
(.venv) (base) [krenzm@justhpc-login Python_skripts]$ emacs extract_isosurface_height_map_parallel.py
(.venv) (base) [krenzm@justhpc-login Python_skripts]$ python extract_isosurface_height_map_parallel.py -h
usage: extract_isosurface_height_map_parallel.py [-h] [--scan_start SCAN_START] [--pp_deflection_dir PP_DEFLECTION_DIR] [--offset OFFSET] [--tip_distance TIP_DISTANCE] [--b_decay B_DECAY] [--no_tip_flip_cutoff] [--Amp AMP] [--iso_range ISO_RANGE] [--offset_range OFFSET_RANGE]
                                                 [--decay_range DECAY_RANGE]
                                                 infile_chg outfile iso_value

Extracts the height profile of an electron density isosurface. Includes the probe-particle deflections if provided

positional arguments:
  infile_chg            Charge density in the xsf format
  outfile               Name of the height profile outfile
  iso_value             Iso value at which the height will be extracted

optional arguments:
  -h, --help            show this help message and exit
  --scan_start SCAN_START, -s SCAN_START
                        relative height (from 0 to 1) where the height scan begins. Sometimes necessary if there are artifacts in the density
  --pp_deflection_dir PP_DEFLECTION_DIR, -ppd PP_DEFLECTION_DIR
                        Directory where the probe particle deflections are stored. Deflections have to have the filenames PPpos_x.xsf, PPpos_y.xsf and PPpos_z.xsf, as in the ppafm default. If no directory is given, deflections are not included (pure isosurface scan)
  --offset OFFSET, -o OFFSET
                        z-offset between the probe particle deflections and the charge density. If the constant height AFM images are adjusted to fit the charge density, this is not necessary. Can be useful to find the correct parameters for the AFM calculation.
  --tip_distance TIP_DISTANCE, -R TIP_DISTANCE
                        Equilibrium distance between probe-particle and tip. To be chosen as in ppafm
  --b_decay B_DECAY, -b B_DECAY
                        Decay constant for tunneling current changes due to tip-pp distance changes (shorter distance more current, longer distance less current).
  --no_tip_flip_cutoff, -tfc
                        Disables the detection of unphysical probe particle flips (probe particle higher than tip). This should not be necessary, used only for errorshooting
  --Amp AMP, -A AMP     Peak-to-Peak amplitude over which the charge density is averaged. Should be the same as in ppafm
  --iso_range ISO_RANGE, -ir ISO_RANGE
                        range of iso values to be used. Given either with a , separated list (e.g 0.00002, 0.00003, etc.) or as a range with syntax "start stop stepsize". Overwrites iso_value.
  --offset_range OFFSET_RANGE, -or OFFSET_RANGE
                        range of offsets to be used. Given either with a , separated list (e.g 0.4, 0.5, etc.) or as a range with syntax "start stop stepsize". Overwrites offset.
  --decay_range DECAY_RANGE, -dr DECAY_RANGE


                        range of tip-pp decays to be used. Given either with a , separated list (e.g 0.05, 0.06, etc.) or as a range with syntax "start stop stepsize". Overwrites b_decay.
