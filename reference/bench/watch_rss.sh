#!/bin/bash
# Sample RSS of the segmenter every 2s until it exits.
"$@" > run15k.log 2>&1 &
pid=$!
echo "pid=$pid" > rss15k.txt
peak=0
while kill -0 $pid 2>/dev/null; do
    rss=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' ')
    [ -z "$rss" ] && break
    mb=$((rss / 1024))
    [ "$mb" -gt "$peak" ] && peak=$mb
    echo "$(date +%s) ${mb}MB" >> rss15k.txt
    sleep 2
done
wait $pid; rc=$?
echo "exit=$rc peak=${peak}MB" >> rss15k.txt
echo "exit=$rc peak=${peak}MB"
