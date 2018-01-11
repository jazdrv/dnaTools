#!/bin/sh

cd $REDUX_ENV
sqlite3 variant.db<<!
.fullschema
!
