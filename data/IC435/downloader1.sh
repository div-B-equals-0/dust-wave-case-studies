#! /bin/sh
echo  '************************************************'
echo  '*'
echo  '*    Date: Fri Jul 26 10:37:44 PDT 2019'
echo  '*'
echo  '*    Download Spitzer data from IRSA'
echo  '*'
echo  '************************************************'

echo;echo  '>> downloading HD38087-6-selected_SMs.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid64162628112_shaapp2_0.zip&return=HD38087-6-selected_SMs.zip&log=true&track=true" --output-document=HD38087-6-selected_SMs.zip
echo;echo  '>> downloading HD38087-6-selected_SMs-part2.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid64162628112_shaapp2_1.zip&return=HD38087-6-selected_SMs-part2.zip&log=true&track=true" --output-document=HD38087-6-selected_SMs-part2.zip
echo; echo; echo '*** All downloads and extractions (if requested) completed ***'
