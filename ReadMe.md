# Usage Instructions

Users can utilize this script by modifying the start date and the duration of the download within the main function of each script.

## Directory Structure

- **NavDecoder/bin**: Contains the necessary binary executables for downloading GNSS products.
- **NavDecoder/download**: Contains the required scripts for downloading various GNSS data.

## Supported Products

1. **Broadcast Navigation Files**: Supports Rinex 3.x and Rinex 4.x formats.
2. **Precise Orbits/Clock Files**:
   - **Near-Real-Time (NRT)**:
     - B2B: [PPP-B2b]
     - HAS: [Galileo-HAS]
     - WUM_NRT: [WUM0MGXNRT_]
   - **Rapid**:
     - CNT: [CNES Real-time archive]
     - GFR: [GFZ0MGXRAP_]
     - WHR: [WUM0MGEXRAP_]
   - **Final**:
     - GRM: [GRG0MGEXFINX_]
     - IGS: [IGS combined solution]
     - WUM: [WUM0MGEXFINX_]
3. **DCB Products**: Provided by CAS.
4. **IGS Sinex Coordinates File**: Contains station coordinate information.
5. **Rinex Navigation Files**: Supports both Rinex 3.x and 4.x formats.

## Python Script Functions

- **cmn_tools.py**: Offers essential utilities for time and coordinate conversion.
- **down_eph_clk.py**: Downloads orbit and clock data from IGS/MGEX data centers.
- **down_PPP_products.py**: Downloads navigation Rinex files (3.x and 4.x), orbits/clocks, and ISB from data centers.
- **down_rinex_DSB.py**: Downloads DCB products from CAS.
- **down_rinex_nav_3x.py**: Downloads Rinex navigation files in 3.x format.
- **down_rinex_nav_4x.py**: Downloads Rinex navigation files in 4.x format.
- **down_rinex_obs.py**: Downloads Rinex observation files.
- **down_sinex_pos.py**: Downloads Sinex Pos files and extracts reference coordinates for the specified station.

## Usage Instructions

To begin using the scripts, users should modify the start date and the length of the download period by adjusting the parameters in the main function of each respective script.
