#!/bin/sh

g++ -I./ -isystem eigen_3.2.2 -isystem boost_1.54.0 -Wall -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -pipe -Wno-unused-function -ftemplate-depth-256   -c -O3 -o var_stack.o stan/agrad/rev/var_stack.cpp
ar -rs libstan.a  var_stack.o
