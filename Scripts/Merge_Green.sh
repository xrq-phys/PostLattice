#!/bin/bash

if [ -z $3 ]; then
	echo "Usage: Merge_Green.sh Green1.def Green2.def GreenOut.def"
fi

linetot=$(expr $(cat $1 | sed '/^\s*$/d' | wc -l) + $(cat $2 | sed '/^\s*$/d' | wc -l) - 5)

cat > $3 << EOF
-------------------
 GreenOut $linetot
-------------------
-------------------
-------------------
EOF

cat $1 | sed 1d | sed 1d | sed 1d | sed 1d | sed 1d >> $3
cat $2 | sed 1d | sed 1d | sed 1d | sed 1d | sed 1d >> $3

echo "Done. Please merge the output files then."

