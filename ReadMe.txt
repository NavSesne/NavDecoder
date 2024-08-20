Usage:
  User can use the script by change the strat data, length of the downloading from the main of each script

NavDecoder/bin: Contains the necessary binary executables for downloading GNSS products.
NavDecoder/download: Houses the required scripts for downloading various GNSS data.

Python Script Functions
cmn_tools.py: Provides essential tools for time and coordinate conversion.
down_eph_clk.py: Downloads orbits and clocks from the IGS/MGEX data centers.
down_PPP_products.py: Downloads navigation Rinex 3.X and 4.X files, orbits/clocks, and ISB from data centers.
down_rinex_DSB.py: Downloads the DCB products from CAS.
down_rinex_nav_3x.py: Downloads the Rinex navigation 3.X file.
down_rinex_nav_4x.py: Downloads the Rinex navigation 4.X file.
down_rinex_obs.py: Downloads the Rinex observation file.
down_sinex_pos.py: Downloads the Sinex Pos file and retrieves the reference coordinates for the specified station.
Usage Instructions
Users can modify the start date and the duration of the download by adjusting the parameters in the main function of each script.

