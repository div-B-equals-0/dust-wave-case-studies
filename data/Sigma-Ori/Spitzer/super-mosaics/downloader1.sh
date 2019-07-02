#! /bin/sh
echo  '************************************************'
echo  '*'
echo  '*    Date: Mon Jul 01 21:25:05 PDT 2019'
echo  '*'
echo  '*    Download Spitzer data from IRSA'
echo  '*'
echo  '************************************************'

echo;echo  '>> downloading sigmaOri-20-selected_SMs.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid62041352915_shaapp1_0.zip&return=sigmaOri-20-selected_SMs.zip&log=true&track=true" --output-document=sigmaOri-20-selected_SMs.zip
echo;echo  '>> downloading sigmaOri-20-selected_SMs-part2.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid62041352915_shaapp1_1.zip&return=sigmaOri-20-selected_SMs-part2.zip&log=true&track=true" --output-document=sigmaOri-20-selected_SMs-part2.zip
echo;echo  '>> downloading sigmaOri-20-selected_SMs-part3.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid62041352915_shaapp1_2.zip&return=sigmaOri-20-selected_SMs-part3.zip&log=true&track=true" --output-document=sigmaOri-20-selected_SMs-part3.zip
echo;echo  '>> downloading sigmaOri-20-selected_SMs-part4.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid62041352915_shaapp1_3.zip&return=sigmaOri-20-selected_SMs-part4.zip&log=true&track=true" --output-document=sigmaOri-20-selected_SMs-part4.zip
echo; echo; echo '*** All downloads and extractions (if requested) completed ***'
