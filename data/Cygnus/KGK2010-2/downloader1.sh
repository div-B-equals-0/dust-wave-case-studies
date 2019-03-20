#! /bin/sh
echo  '************************************************'
echo  '*'
echo  '*    Date: Tue Mar 19 08:03:22 PDT 2019'
echo  '*'
echo  '*    Download Spitzer data from IRSA'
echo  '*'
echo  '************************************************'

echo;echo  '>> downloading 2MASSJ20343455%2B4158295-18-selected_SMs.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid53007615793_shaapp1_0.zip&return=2MASSJ20343455%2B4158295-18-selected_SMs.zip&log=true&track=true" --output-document=2MASSJ20343455%2B4158295-18-selected_SMs.zip
echo;echo  '>> downloading 2MASSJ20343455%2B4158295-18-selected_SMs-part2.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid53007615793_shaapp1_1.zip&return=2MASSJ20343455%2B4158295-18-selected_SMs-part2.zip&log=true&track=true" --output-document=2MASSJ20343455%2B4158295-18-selected_SMs-part2.zip
echo;echo  '>> downloading 2MASSJ20343455%2B4158295-18-selected_SMs-part3.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid53007615793_shaapp1_2.zip&return=2MASSJ20343455%2B4158295-18-selected_SMs-part3.zip&log=true&track=true" --output-document=2MASSJ20343455%2B4158295-18-selected_SMs-part3.zip
echo;echo  '>> downloading 2MASSJ20343455%2B4158295-18-selected_SMs-part4.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid53007615793_shaapp1_3.zip&return=2MASSJ20343455%2B4158295-18-selected_SMs-part4.zip&log=true&track=true" --output-document=2MASSJ20343455%2B4158295-18-selected_SMs-part4.zip
echo;echo  '>> downloading 2MASSJ20343455%2B4158295-18-selected_SMs-part5.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid53007615793_shaapp1_4.zip&return=2MASSJ20343455%2B4158295-18-selected_SMs-part5.zip&log=true&track=true" --output-document=2MASSJ20343455%2B4158295-18-selected_SMs-part5.zip
echo;echo  '>> downloading 2MASSJ20343455%2B4158295-18-selected_SMs-part6.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid53007615793_shaapp1_5.zip&return=2MASSJ20343455%2B4158295-18-selected_SMs-part6.zip&log=true&track=true" --output-document=2MASSJ20343455%2B4158295-18-selected_SMs-part6.zip
echo;echo  '>> downloading 2MASSJ20343455%2B4158295-18-selected_SMs-part7.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid53007615793_shaapp1_6.zip&return=2MASSJ20343455%2B4158295-18-selected_SMs-part7.zip&log=true&track=true" --output-document=2MASSJ20343455%2B4158295-18-selected_SMs-part7.zip
echo;echo  '>> downloading 2MASSJ20343455%2B4158295-18-selected_SMs-part8.zip ...'
wget "https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid53007615793_shaapp1_7.zip&return=2MASSJ20343455%2B4158295-18-selected_SMs-part8.zip&log=true&track=true" --output-document=2MASSJ20343455%2B4158295-18-selected_SMs-part8.zip
echo; echo; echo '*** All downloads and extractions (if requested) completed ***'
