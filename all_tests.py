import sys

sys.path.insert(0, './tests')

print('Running general test...')
import test
test.runTest()
print('Done!')


print('Running parsing test...')
import parsing_test
parsing_test.runTest()
print('Done!')


print('Running propagation test...')
import propagation_test
propagation_test.runTest()
print('Done!')
