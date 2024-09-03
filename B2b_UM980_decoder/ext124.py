# Let's extract all content that starts with #PPPB2BINFO from the file

def extract_data_from_line1(line):
    # 将每一行按逗号分割
    parts = line.split(',')

    # 检查是否有足够的数据
    if len(parts) < 15:
        return None

    # 提取所需的数据
    type = parts[0]
    week = parts[4]
    wn = str(int(parts[5]) // 1000)
    prn = str(int(parts[10]) - 160)
    iodssr = parts[11]
    iodp = parts[12]
    todbdt = parts[13]

    # 将第15个数据转换为二进制，并去除'*'后的部分，确保长度为原十六进制字符串长度的4倍
    hex_data = parts[14].split('*')[0]
    binary_data = bin(int(hex_data, 16))[2:].zfill(len(hex_data) * 4)

    # 切片提取不同卫星系统的数据
    bds = binary_data[:63]
    gps = binary_data[63:100]
    galileo = binary_data[100:137]
    glonass = binary_data[137:174]

    return f'{type} week:{week} wn:{wn} prn:{prn} iodssr:{iodssr} iodp:{iodp} todbdt:{todbdt} BDS:{bds} GPS:{gps} Galileo:{galileo} GLONASS:{glonass}\n'

def extract_data_from_line2(line):
    parts = line.split(',')
    if len(parts) < 56:  # 确保每行有足够的数据
        return None
    type = parts[0]
    week=parts[4]
    wn = str(int(parts[5]) // 1000)
    prn = str(int(parts[10]) - 160)
    iodssr = parts[11]
    iodp = parts[12]
    todbdt = parts[13]
    result = f'{type} week:{week} wn:{wn} prn:{prn} iodssr:{iodssr} iodp:{iodp} todbdt:{todbdt}'

    for i in range(14, len(parts) - 1, 7):
        satslot = parts[i]
        iodn = parts[i + 1]
        Rorb = float(parts[i + 2]) * 0.0016
        Aorb = float(parts[i + 3]) * 0.0064
        Corb = float(parts[i + 4]) * 0.0064
        iodcorr = parts[i + 5]
        URAI = parts[i + 6].split('*')[0]  # 去除可能存在的'*'符号
        result += f' satslot:{satslot} iodn:{iodn} Rorb:{Rorb:.4f} Aorb:{Aorb:.4f} Corb:{Corb:.4f} iodcorr:{iodcorr} URAI:{URAI}'

    return result + '\n'

def extract_data_from_line4(line):
    parts = line.split(',')
    if len(parts) < 56:  # 确保每行有足够的数据
        return None
    type = parts[0]
    week = parts[4]
    wn = str(int(parts[5]) // 1000)
    prn = str(int(parts[10]) - 160)
    iodssr = parts[11]
    iodp = parts[12]
    todbdt = parts[13]
    sub = parts[14]
    result = f'{type} week:{week} wn:{wn} prn:{prn} sub:{sub} iodssr:{iodssr} iodp:{iodp} todbdt:{todbdt}'
    for i in range(18, len(parts) - 1, 2):
        iodcorr = parts[i]
        sc0 = float(parts[i + 1]) * 0.0016
        result += f' iodcorr:{iodcorr} sc0:{sc0}'

    return result + '\n'
def extract_all_pppb2binfo_content(file_path, output_file_path):
    with open(file_path, 'rb') as f:
        lines = f.read()
        parts = lines.split(b'#PPPB2BINFO')[1:]  # Skip the first part as it doesn't start with #PPPB2BINFO
        with open(output_file_path, 'w') as output_file:
            for part in parts:
                # Extract content up to the next newline or the end of the file
                end_index = part.find(b'\n')
                if end_index != -1:
                    content = part[:end_index]
                else:
                    content = part

                # Decode, replace semicolons with commas, and remove 'A' from the first value
                decoded_line = content.decode('ascii', errors='ignore')
                replaced_line = decoded_line.replace(';', ',').replace('*', ',')
                first_value, rest = replaced_line.split(',', 1)
                first_value = first_value.replace('A', '')  # Remove 'A' from the first value
                processed_line = first_value + ',' + rest

                # Call the appropriate function based on the first value
                if first_value == '1':
                    extracted_data = extract_data_from_line1(processed_line)
                elif first_value == '2':
                    extracted_data = extract_data_from_line2(processed_line)
                elif first_value == '4':
                    extracted_data = extract_data_from_line4(processed_line)
                else:
                    extracted_data = None

                # Write the extracted data to the output file
                if extracted_data:
                    output_file.write(extracted_data)


# # Call the function with the input and output file paths
# input_file_path = 'D:\\cssrlib-main\\src\\cssrlib\\um982_raw_data\\log_UM982_20240319_00.txt'
# output_file_path = 'D:\\cssrlib-main\\src\\cssrlib\\um982_raw_data\\0319.txt'
# extract_all_pppb2binfo_content(input_file_path, output_file_path)





