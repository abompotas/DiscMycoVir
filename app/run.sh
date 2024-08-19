#!/bin/bash

cp -r /virus-discovery/* /usr/share/nginx/html
sed -i 's@___discvirAPI___@'"${DISCVIR_API}"'@g' /usr/share/nginx/html/main.js
nginx -g 'daemon off;'
