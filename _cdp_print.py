import base64, json, subprocess, time, pathlib
import requests
import websocket

here = pathlib.Path(__file__).parent
html_path = here / "design_report.html"
pdf_path = here / "design_report.pdf"

port = 9444
proc = subprocess.Popen([
    "chromium-browser", "--headless=new", "--disable-gpu", "--no-sandbox",
    "--user-data-dir=/tmp/chrome-profile-cdp",
    f"--remote-debugging-port={port}", "--remote-allow-origins=*",
    "about:blank",
])
time.sleep(2)

tabs = None
for _ in range(30):
    try:
        tabs = requests.get(f"http://127.0.0.1:{port}/json").json()
        if tabs:
            break
    except Exception:
        pass
    time.sleep(0.5)

pages = [t for t in tabs if t.get("type") == "page"]
tab = pages[0] if pages else tabs[0]
ws_url = tab["webSocketDebuggerUrl"]
ws = websocket.create_connection(ws_url, timeout=30)

def send(method, params=None, msg_id=1):
    ws.send(json.dumps({"id": msg_id, "method": method, "params": params or {}}))
    while True:
        resp = json.loads(ws.recv())
        if resp.get("id") == msg_id:
            return resp

send("Page.enable", msg_id=1)
nav = send("Page.navigate", {"url": f"file://{html_path}"}, msg_id=2)
print("navigate:", nav)

time.sleep(6)  # let MathJax typeset

result = send("Page.printToPDF", {
    "printBackground": True,
    "displayHeaderFooter": False,
    "preferCSSPageSize": True,
    "marginTop": 0.5, "marginBottom": 0.5, "marginLeft": 0.5, "marginRight": 0.5,
}, msg_id=3)

if "result" not in result:
    print("ERROR:", result)
else:
    data = result["result"]["data"]
    pdf_path.write_bytes(base64.b64decode(data))
    print(f"wrote {pdf_path} ({pdf_path.stat().st_size} bytes)")

ws.close()
proc.terminate()
