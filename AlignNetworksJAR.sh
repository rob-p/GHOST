#!/usr/bin/env bash

SYSTEM=$(uname)

if [ "$SYSTEM" = "Linux" ] 
then
    MEM=$(cat /proc/meminfo | grep MemTotal | awk '{print $2}')
    MEM=$(echo "scale=2; ${MEM}/2^20" | bc)
elif [ "$SYSTEM" = "Darwin" ]
then
    MEM=$(sysctl hw.memsize | awk '{print $2}')
    MEM=$(echo "scale=2; ${MEM}/2^30" | bc)
else
    echo "OSX and Linux are the only supported platforms"
fi

echo "Running on ${SYSTEM} with $MEM GB of total RAM"

JVM_MAX=$(printf "%.0f" $(echo "${MEM}*0.95" | bc))

ARGS=$@

JAVA_OPTS="-Xmx${JVM_MAX}G -XX:+UseConcMarkSweepGC -XX:+CMSClassUnloadingEnabled" 

java ${JAVA_OPTS} -jar Ghost.jar ${ARGS}