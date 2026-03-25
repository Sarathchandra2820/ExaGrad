name=$(git remote)
echo $name

git add .
git commit -m "auto changes"
git push -u $name main

