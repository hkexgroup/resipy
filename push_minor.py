# run with python push_minor.py

import re
import subprocess
import os

# get branch name
branch = subprocess.check_output(['git', 'branch', '--show-current']).decode().strip()

# merge master to avoid merge conflicts later
print('===== Fetch master, then merge it ====')
os.system('git fetch origin master')
os.system('git merge origin/master')

print("=== Update version number ===")
# increment version number
with open('src/resipy/Project.py', 'r') as f:
    x = f.read()

out = re.findall("(\d+\.\d+\.\d+)", x)
print(out)
version_old = out[0]
version = version_old.split('.')
version[-1] = str(int(version[-1])+1)
version_new = '.'.join(version)

with open('src/resipy/Project.py', 'w') as f:
    f.write(x.replace(version_old, version_new))
   
with open('src/setup.py', 'rw') as f:
    f.write(f.read().replace(version_old, version_new))

with open('src/Version.details', 'rw') as f:
    f.write(f.read().replace(version_old, version_new).replace(','.join(version_old.split('.')), ','.join(version_new.split('.'))))

# commit version increase
os.system('git add src/resipy/Project.py')
os.system(f'git commit -m "Update to v{version_new}"')

print("=== Create tag and merge request ===")
# create tag
os.system(f'git tag v{version_new}')

# push and create merge request
os.system(f'git push --set-upstream origin {branch} -o merge_request.create -o merge_request.title="v{version_new}" tag v{version_new}')


