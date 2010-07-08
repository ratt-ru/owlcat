#!/bin/bash
# -*- coding: utf-8 -*-
if [ "$2" == "" ]; then
  echo "Usage: $0 msname ifr_subset_string"
  exit 2
fi

python -c "import Owlcat;import Meow.IfrSet;ss=Meow.IfrSet.from_ms('$1').subset('$2'); print '%d//%s'%(len(ss.ifr_index()),ss.taql_string());" 2>/dev/null | grep ANTENNA1==
