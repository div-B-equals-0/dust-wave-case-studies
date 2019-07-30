#! /bin/sh
echo  '************************************************'
echo  '*'
echo  '*    Date: Thu Jul 25 11:01:16 PDT 2019'
echo  '*'
echo  '*    Download Spitzer data from IRSA'
echo  '*'
echo  '************************************************'

echo;echo  '>> downloading sigmaori-19-selected_SMs.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid64077572981_shaapp2_0.zip&return=sigmaori-19-selected_SMs.zip&log=true&track=true" --output-document=sigmaori-19-selected_SMs.zip
echo;echo  '>> downloading sigmaori-19-selected_SMs-part2.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid64077572981_shaapp2_1.zip&return=sigmaori-19-selected_SMs-part2.zip&log=true&track=true" --output-document=sigmaori-19-selected_SMs-part2.zip
echo;echo  '>> downloading sigmaori-19-selected_SMs-part3.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid64077572981_shaapp2_2.zip&return=sigmaori-19-selected_SMs-part3.zip&log=true&track=true" --output-document=sigmaori-19-selected_SMs-part3.zip
echo; echo; echo '*** All downloads and extractions (if requested) completed ***'
