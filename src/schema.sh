#!/bin/bash

cd $REDUX_ENV

DB="variant"
TB="$1"

echo ""
sqlite3 $DB.db <<!
.schema $TB
!
echo ""
