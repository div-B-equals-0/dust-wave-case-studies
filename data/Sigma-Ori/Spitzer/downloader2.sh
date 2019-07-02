#! /bin/sh
echo  '************************************************'
echo  '*'
echo  '*    Date: Mon Jul 01 21:38:27 PDT 2019'
echo  '*'
echo  '*    Download Spitzer data from IRSA'
echo  '*'
echo  '************************************************'

echo;echo  '>> downloading sigmaOri-2-selected_Post_BCDs.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid62042220674_shaapp1_0.zip&return=sigmaOri-2-selected_Post_BCDs.zip&log=true&track=true" --output-document=sigmaOri-2-selected_Post_BCDs.zip
echo;echo  '>> downloading sigmaOri-2-selected_Post_BCDs-part2.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid62042220674_shaapp1_1.zip&return=sigmaOri-2-selected_Post_BCDs-part2.zip&log=true&track=true" --output-document=sigmaOri-2-selected_Post_BCDs-part2.zip
echo;echo  '>> downloading sigmaOri-2-selected_Post_BCDs-part3.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid62042220674_shaapp1_2.zip&return=sigmaOri-2-selected_Post_BCDs-part3.zip&log=true&track=true" --output-document=sigmaOri-2-selected_Post_BCDs-part3.zip
echo; echo; echo '*** All downloads and extractions (if requested) completed ***'
