cat result_fixed.html | \
sed -e 's/tara_oceans/Tara_oceans/g' | \
sed -e 's/muddy/Muddy/g' | \
sed -e 's/sandy/Sandy/g' | \
sed -e 's/moose/Moose/g' | \
sed -e 's/atlantic/Atlantic/g' | \
sed -e 's/plant_microflora/Plant_microflora/g' | \
sed -e 's/sludge_soil/Sludge_soil/g' | \
sed -e 's/cow_dung/Cow_dung/g' | \
sed -e 's/mangrove/Mangrove/g' | \
sed -e 's/gose/GOS/g' | \
sed -e 's/human_faeces/Human_faeces/g' | \
sed -e 's/skin/Skin/g' > result_fixed2.html
