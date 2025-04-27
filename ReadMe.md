# Usage Instructions

This Python toolbox is designed for **BDS PPP-B2b** and **Galileo HAS decoding**. The key features of the **NavDecoder** are as follows:

1. **Download tools** for retrieving navigation and MGEX orbit/clock products.
2. **Support for raw binary HAS data** from Septentrio GNSS receivers. The system can be easily extended to support binary data from various manufacturers by extracting raw Galileo C/NAV data for input into the decoding function.
3. **Support for B2b data** from both **Septentrio** and **Unicore** GNSS receivers:
   - Septentrio provides raw binary format, which requires an LDPC decoder to recover plain corrections.
   - Unicore offers a plain ASCII format that can be directly processed to recover corrections.
4. **Capability to save corrections** in the **BNC universal format**.
5. **Capability to save corrections** in **SP3** and **CLK formats**.
6. **Provision of an archive** for **B2b** and **HAS corrections**.

### Main Functions

- **Decoding the Septentrio HAS Binary format** from the `"GALRawCNAV"` block: `decode_HAS_sept.py`
- **Decoding the Septentrio B2b Binary format** from the `"BDSRawB2b"` block: `decode_B2B_sept.py`
- **Decoding the Unicor B2b Binary format** from the `"PPPB2bInfo"` block: `decode_B2B_UM982.py`

### Configuration for the decoding:

The configuration settings are located in the following section of the code:

```python
# Start of configuration
max_orbit_delay = 300  # Maximum orbit delay in seconds
max_clock_delay = 30   # Maximum clock delay in seconds
start_date = datetime(2024, 5, 14)  # Start date for the process
process_days = 1  # Number of days to process
file_has_template = r'test_data\SEPT{}0.{}__SBF_GALRawCNAV.txt'  # Template for HAS files
nav_file_template = r'test_data\eph\BRDC00IGS_R_{}0000_01D_MN.rnx'  # Template for navigation files
corr_dir_template = r'test_data\SEPT{}_HAS'  # Directory for correction files
# End of configuration
 ```

You only need to modify the **start date**, **processing duration**, and **file paths**.
## Directory Structure
- **NavDecoder/B2b_HAS_decoder**: Tools for decoding Septentrio B2b HAS data.
- **NavDecoder/B2b_UM980_decoder**: Tools for decoding Unicore B2b data.
- **NavDecoder/bin**: Contains binary executables for downloading GNSS products.
- **NavDecoder/download**: Contains scripts for downloading various GNSS data.
- **NavDecoder/test_data**: Example data for decoding.

## Supported Products of the download tool
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


# PPP-B2b and Galileo HAS SSR Correction Data Archive

This repository provides download links for the PPP-B2b and Galileo HAS SSR correction data archives.

## Download Links

### **navcm_HAS_ssr**

- **Download Link:** [Baidu Netdisk](https://pan.baidu.com/s/1yv14UxF68HHLWYyZXJ3aEA)  
- **Extraction Code:** `7wn6`

---

### **navcm_B2b_ssr**

- **Download Link:** [Baidu Netdisk](https://pan.baidu.com/s/1DdUYe-AnFMM-5id4bkjt0Q)  
- **Extraction Code:** `6ypj`

---

### **SEPT RAW files**

- **Download Link:** [Baidu Netdisk]( https://pan.baidu.com/s/1d-a8jYOrjOwgWNmt5PqWmQ)  
- **Extraction Code:** `h3cf`

---

### **UM98 RAW files**

- **Download Link:** [Baidu Netdisk](https://pan.baidu.com/s/1zZNejuNJoinoNzsq8-6pMA )  
- **Extraction Code:** `nft7`

---
## How to Access

1. Click on the provided Baidu Netdisk download links.
2. Enter the corresponding extraction code when prompted.
3. Download and extract the archive files.
