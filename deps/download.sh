#!/bin/bash

odlcommontools_git="https://github.com/OptimalDesignLab/ODLCommonTools.jl.git"
if [ ! -d "ODLCommonTools" ]; then
  git clone $odlcommontools_git ../../ODLCommonTools
fi


