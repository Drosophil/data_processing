'''
It is not a part of a library, it's an utility
for making from HUGE_Data_Set.csv small sample.csv samples
'''

with open('Practice/Practice/HUGE_Data_Set.csv', 'r') as f:
    with open('sample.csv', 'w') as g:
        for i in range(30000):
            g.write(f.readline())
