#!/bin/bash

sleep 1
echo 'hi from stdout'
sleep 1
echo 'hi from stderr' >&2
