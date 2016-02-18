#!/bin/bash

(make && ./main) || exit 1

(diff channel.dat goldenfile && echo "PASSED") || echo "FAILED: goldenfile has changed!"

