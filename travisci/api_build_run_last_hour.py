import sys
import time
from datetime import datetime, date
import json

def api_build_run_last_hour(x):
    return x['state'] != 'canceled' and \
                         x['event_type'] == 'api' and \
                         (x['finished_at'] == None or \
                          float((datetime.fromtimestamp(time.time()) - datetime.strptime(x['finished_at'],'%Y-%m-%dT%H:%M:%SZ')).total_seconds())/3600 < 1.0)

builds = list(filter(api_build_run_last_hour, json.load(sys.stdin)['builds']))
print(len(builds) > 0)
