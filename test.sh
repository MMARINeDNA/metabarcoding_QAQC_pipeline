#!/bin/sh


while true; do
read -p 'Another round of Taxonomy? (y/n) ' tax_reply
if [[ ${tax_reply} == "y" ]]; then
cat assignTaxonomy_codes
read -p "Enter Primer Name: " primer

# mifish
if [[ $primer = "MFU" ]]
then
echo mifish is done!

# marver 1
elif [[ $primer = "MV1" ]]
then
echo mv1 is done!


# marver 3
elif [[ $primer = "MV3" ]]
then
echo mv3 is done!


# d-loop
elif [[ $primer = "DL" ]]
then
echo d loop is done!


# leray
elif [[ $primer = "LRY" ]]
then
echo leray is done!


# C16
elif [[ $primer = "C16" ]]
then
echo c16 is done!


# C18 
elif [[ $primer = "C18" ]]
then
echo c18 is done!

else 
echo idk that one... did you make a typo?
fi

else
break
fi
done