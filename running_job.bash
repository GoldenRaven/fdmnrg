#!/bin/bash
expect running_job.expect |grep -v spawn|grep -v login > /tmp/tmpp
cat /tmp/tmpp
