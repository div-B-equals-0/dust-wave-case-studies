#! /bin/sh
echo  '************************************************'
echo  '*'
echo  '*    Date: Tue Mar 19 08:45:26 PDT 2019'
echo  '*'
echo  '*    Download Spitzer data from IRSA'
echo  '*'
echo  '************************************************'

echo;echo  '>> downloading 2MASSJ20343455%2B4158295-5-selected_SMs.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid53010149204_shaapp1_0.zip&return=2MASSJ20343455%2B4158295-5-selected_SMs.zip&log=true&track=true" --output-document=2MASSJ20343455%2B4158295-5-selected_SMs.zip
echo;echo  '>> downloading 2MASSJ20343455%2B4158295-5-selected_SMs-part2.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid53010149204_shaapp1_1.zip&return=2MASSJ20343455%2B4158295-5-selected_SMs-part2.zip&log=true&track=true" --output-document=2MASSJ20343455%2B4158295-5-selected_SMs-part2.zip
echo; echo; echo '*** All downloads and extractions (if requested) completed ***'
