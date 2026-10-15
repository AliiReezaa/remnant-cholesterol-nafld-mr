#!/bin/bash
# usage: cr.sh queries_file ; prints top-1 Crossref hit per line
while IFS= read -r t; do
  [ -z "$t" ] && continue
  echo "## $t" | cut -c1-120
  curl -s -m 25 -G "https://api.crossref.org/works" --data-urlencode "query.bibliographic=$t" --data-urlencode "rows=1" --data-urlencode "select=DOI,title,author,container-title,issued,volume,issue,page,article-number,type" | $HOME/mrenv/bin/python -I -c "
import sys,json
try:
  for m in json.load(sys.stdin)['message']['items']:
    au=m.get('author') or []; y=(m.get('issued',{}).get('date-parts') or [[None]])[0][0]
    names=', '.join((a.get('family','')+' '+''.join(w[0] for w in a.get('given','').replace('-',' ').split())) for a in au[:3])+(' et al' if len(au)>3 else '')
    print('  ',m['DOI'],'|',names,'|',(m.get('title') or [''])[0][:110],'|',(m.get('container-title') or [''])[0],'|',y,'|',m.get('volume'),m.get('issue'),m.get('page') or m.get('article-number'),'|',m.get('type'))
except Exception as e: print('  FAIL',e)"
  sleep 0.4
done < "$1"
