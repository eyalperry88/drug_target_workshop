import os, sys

print 'Running general test...'
from tests import test
print 'Done!'


print 'Running parsing test...'
from tests import parsing_test
print 'Done!'


print 'Running propagation test...'
from tests import propagation_test
print 'Done!'
