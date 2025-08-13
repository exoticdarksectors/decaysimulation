import re

with open('dy_cross.txt', 'w') as output_file:
    # Open the HTML file and read line by line
    with open('DY-output-test1/crossx.html', 'r') as file:
        for line in file:
            # Find and extract characters between 'results.html">' and '<font face=symbol>&#177'
            pattern = r'results.html">(.+?)<font face=symbol>&#177'
            match = re.search(pattern, line)
            if match:
                extracted_text = match.group(1)  # Extract the matched characters
                # Convert the extracted text to double-type numbers
                try:
                    extracted_number = float(extracted_text)
                    # Write the number to the output file
                    output_file.write(f'{extracted_number}\n')
                except ValueError:
                    pass  # Handle conversion errors
