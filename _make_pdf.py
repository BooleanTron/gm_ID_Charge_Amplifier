import markdown, pathlib

here = pathlib.Path(__file__).parent
md_text = (here / "design_report.md").read_text()

body = markdown.markdown(md_text, extensions=["tables", "fenced_code", "sane_lists"])

html = f"""<!doctype html>
<html><head><meta charset="utf-8">
<script>
window.MathJax = {{
  tex: {{ inlineMath: [['$', '$']], displayMath: [['$$', '$$']] }},
  svg: {{ fontCache: 'global' }}
}};
</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>
<style>
  @page {{ size: A4; margin: 20mm 18mm; }}
  body {{ font-family: Georgia, 'Times New Roman', serif; font-size: 10.5pt;
         line-height: 1.5; color: #1a1a1a; max-width: 100%; }}
  h1 {{ font-size: 18pt; border-bottom: 2px solid #333; padding-bottom: 6px; }}
  h2 {{ font-size: 14pt; margin-top: 1.6em; border-bottom: 1px solid #999; padding-bottom: 3px; }}
  h3 {{ font-size: 12pt; margin-top: 1.3em; color: #223; }}
  h4 {{ font-size: 11pt; margin-top: 1em; }}
  table {{ border-collapse: collapse; width: 100%; margin: 0.8em 0; font-size: 9.3pt; }}
  th, td {{ border: 1px solid #888; padding: 4px 7px; text-align: left; }}
  th {{ background: #eef0f3; }}
  code {{ background: #f2f2f2; padding: 1px 4px; border-radius: 3px; font-size: 9.5pt; }}
  pre {{ background: #f6f6f6; padding: 8px; border-radius: 4px; overflow-x: auto; }}
  img {{ max-width: 100%; display: block; margin: 10px auto; }}
  blockquote {{ border-left: 3px solid #aaa; margin: 0.5em 0; padding: 2px 12px; color: #444; background:#f9f9f9;}}
  hr {{ border: none; border-top: 1px solid #ccc; margin: 1.4em 0; }}
  strong {{ color: #111; }}
  ul, ol {{ margin: 0.3em 0; padding-left: 1.4em; }}
</style>
</head>
<body>
{body}
</body></html>
"""

(here / "design_report.html").write_text(html)
print("wrote design_report.html")
