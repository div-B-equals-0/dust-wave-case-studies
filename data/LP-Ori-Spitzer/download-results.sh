#! /bin/sh
echo  '************************************************'
echo  '*'
echo  '*    Date: Wed Jun 06 09:47:49 PDT 2018'
echo  '*'
echo  '*    Download Spitzer data from IRSA'
echo  '*'
echo  '************************************************'

echo;echo  '>> downloading LPOri-38-selected_Post_BCDs.zip ...'
curl "http://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid28303618313_shaapp1_0.zip&return=LPOri-38-selected_Post_BCDs.zip&log=true&track=true" -o LPOri-38-selected_Post_BCDs.zip
echo;echo  '>> downloading LPOri-38-selected_Post_BCDs-part2.zip ...'
curl "http://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid28303618313_shaapp1_1.zip&return=LPOri-38-selected_Post_BCDs-part2.zip&log=true&track=true" -o LPOri-38-selected_Post_BCDs-part2.zip
echo;echo  '>> downloading LPOri-38-selected_Post_BCDs-part3.zip ...'
curl "http://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid28303618313_shaapp1_2.zip&return=LPOri-38-selected_Post_BCDs-part3.zip&log=true&track=true" -o LPOri-38-selected_Post_BCDs-part3.zip
echo;echo  '>> downloading LPOri-38-selected_Post_BCDs-part4.zip ...'
curl "http://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid28303618313_shaapp1_3.zip&return=LPOri-38-selected_Post_BCDs-part4.zip&log=true&track=true" -o LPOri-38-selected_Post_BCDs-part4.zip
echo;echo  '>> downloading LPOri-38-selected_Post_BCDs-part5.zip ...'
curl "http://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/Download?file=%24%7Bstage%7D%2Fbid28303618313_shaapp1_4.zip&return=LPOri-38-selected_Post_BCDs-part5.zip&log=true&track=true" -o LPOri-38-selected_Post_BCDs-part5.zip
echo; echo; echo '*** All downloads and extractions (if requested) completed ***'
