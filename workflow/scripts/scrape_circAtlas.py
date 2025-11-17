import requests
import pandas as pd

# The target endpoint
url = "https://ngdc.cncb.ac.cn/circatlas/browse_species1.php"

# Updated headers to match exactly what the browser sends
headers = {
    "User-Agent": "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:134.0) Gecko/20100101 Firefox/134.0",
    "Accept": "application/json, text/javascript, */*; q=0.01",
    "Accept-Language": "en-US,en;q=0.5",
    "Accept-Encoding": "gzip, deflate, br, zstd",
    "Content-Type": "application/x-www-form-urlencoded; charset=UTF-8",
    "X-Requested-With": "XMLHttpRequest",
    "Origin": "https://ngdc.cncb.ac.cn",
    "Connection": "keep-alive",
    "Referer": "https://ngdc.cncb.ac.cn/circatlas/browse1.php"  # Updated referer
}

# Updated payload to match the exact format from the browser request
payload = {
    "draw": "1",
    "columns[0][data]": "species",
    "columns[0][name]": "",
    "columns[0][searchable]": "true",
    "columns[0][orderable]": "true",
    "columns[0][search][value]": "",
    "columns[0][search][regex]": "false",
    "columns[1][data]": "name",
    "columns[1][name]": "",
    "columns[1][searchable]": "true",
    "columns[1][orderable]": "true",
    "columns[1][search][value]": "",
    "columns[1][search][regex]": "false",
    "columns[2][data]": "uid",
    "columns[2][name]": "",
    "columns[2][searchable]": "true",
    "columns[2][orderable]": "true",
    "columns[2][search][value]": "",
    "columns[2][search][regex]": "false",
    "columns[3][data]": "pos",
    "columns[3][name]": "",
    "columns[3][searchable]": "true",
    "columns[3][orderable]": "true",
    "columns[3][search][value]": "",
    "columns[3][search][regex]": "false",
    "columns[4][data]": "strand",
    "columns[4][name]": "",
    "columns[4][searchable]": "true",
    "columns[4][orderable]": "true",
    "columns[4][search][value]": "",
    "columns[4][search][regex]": "false",
    "columns[5][data]": "ctype",
    "columns[5][name]": "",
    "columns[5][searchable]": "true",
    "columns[5][orderable]": "true",
    "columns[5][search][value]": "",
    "columns[5][search][regex]": "false",
    "columns[6][data]": "score",
    "columns[6][name]": "",
    "columns[6][searchable]": "true",
    "columns[6][orderable]": "true",
    "columns[6][search][value]": "",
    "columns[6][search][regex]": "false",
    "columns[7][data]": "nspe",
    "columns[7][name]": "",
    "columns[7][searchable]": "true",
    "columns[7][orderable]": "false",
    "columns[7][search][value]": "",
    "columns[7][search][regex]": "false",
    "columns[8][data]": "ntis",
    "columns[8][name]": "",
    "columns[8][searchable]": "true",
    "columns[8][orderable]": "false",
    "columns[8][search][value]": "",
    "columns[8][search][regex]": "false",
    "columns[9][data]": "nsam",
    "columns[9][name]": "",
    "columns[9][searchable]": "true",
    "columns[9][orderable]": "false",
    "columns[9][search][value]": "",
    "columns[9][search][regex]": "false",
    "columns[10][data]": "len",
    "columns[10][name]": "",
    "columns[10][searchable]": "true",
    "columns[10][orderable]": "true",
    "columns[10][search][value]": "",
    "columns[10][search][regex]": "false",
    "order[0][column]": "0",
    "order[0][dir]": "asc",
    "start": "0",
    "length": "1000",
    "search[value]": "",
    "search[regex]": "false",
    "species": "human"
}

# Storage for all data
all_data = []
total_records = None

# Add session handling
session = requests.Session()

# Make an initial GET request to the main page to get necessary cookies
session.get("https://ngdc.cncb.ac.cn/circatlas/browse1.php")  # Updated initial page

# Set cookies that might be required
cookies = {
    "acceptCookies": "true"
}
session.cookies.update(cookies)

# Add a counter for requests
request_count = 0
max_requests = 780

while True:
    # Send the POST request using the session
    response = session.post(url, headers=headers, data=payload)
    
    # Debug line to check response
    print(f"Request {request_count + 1}: Status Code: {response.status_code}")
    
    response.raise_for_status()
    
    # Parse the JSON response
    json_data = response.json()
    
    # Get total records on first request
    if total_records is None:
        total_records = int(json_data.get("recordsTotal", 0))
        print(f"\nTotal records available: {total_records}")
    
    rows = json_data.get("data", [])
    all_data.extend(rows)
    
    # Calculate progress
    current_count = len(all_data)
    progress = (current_count / total_records) * 100
    print(f"Progress: {current_count}/{total_records} ({progress:.2f}%)")
    
    # Increment counter
    request_count += 1
    
    # Check if we've reached our request limit or if we have all records
    if request_count >= max_requests or current_count >= total_records:
        break

    # Increment for the next page
    payload["start"] = str(int(payload["start"]) + int(payload["length"]))
    payload["draw"] = str(int(payload["draw"]) + 1)

# Print summary
print(f"\nSummary:")
print(f"Total requests made: {request_count}")
print(f"Total entries collected: {len(all_data)}")
print(f"Collection complete: {len(all_data) >= total_records}")

# Convert to a DataFrame and save to CSV
columns = ["species", "name", "uid", "pos", "strand", "ctype", "score", "nspe", "ntis", "nsam", "len"]
df = pd.DataFrame(all_data, columns=columns)
df.to_csv("output.csv", index=False)

print("Data saved to output.csv")
