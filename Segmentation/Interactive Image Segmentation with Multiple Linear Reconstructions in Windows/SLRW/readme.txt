    ======================================================================================================
    Basic Information
    ======================================================================================================
 
    This code is written by Shiming Xiang
 
    Address:    National Laboratory of Pattern Recognition (NLPR), 
                Institute of Automation, Chinese Academy of Sciences, Beijing, 100190, China
    Email:      smxiang@gmail.com

    This software can be used freely for research purposes.
    Published reports of research  using this code (or a modified version) should cite the 
    article that describes the algorithm: 

    Shiming Xiang, Chunhong Pan, Feiping Nie, Changshui Zhang: 
    Interactive Image Segmentation With Multiple Linear Reconstructions in Windows. 
    IEEE Transactions on Multimedia, Vol. 13, No.2, pp. 342-352, 2011.
 

    This code has been compiled and tested using matlab 6.5  and matlab 7.0

    version 1.0, 
    rewritten and updated by Shiming Xiang,  Oct. 19, 2011  


    ======================================================================================================
    WHAT's NEW: !!!!!
    Note that: Our algorithm can also be extended for image matting, although we did not report this in the paper
    In this matlab package, the function for image matting with our method is supported.
    Please run the codes in file: demo_mlrw_interactive_image_matting.m
    ======================================================================================================




    ======================================================================================================
    How to use
    ======================================================================================================

    1. for interactive image segmentation  and interactive image matting, please call function:

       mlrw_for_interactive_image();

       ----------------------------------------------------------------
       The default is designed for  interactive image segmentation. 
       That is, 
       options.type_interaction = 1;
       ----------------------------------------------------------------

    2. There are two demo function. You can directly rum them. Also there are comments in the codes which may be helpful. 
 
       demo_mlrw_interactive_image_segmentation.m
       demo_mlrw_interactive_image_matting.m

      -------------------------------------------------------------------------------------
      Note that, for image matting, we must set: 

      options.type_interaction = 0; 

      and transfer "options" to function "mlrw()", 
      otherwise (default), it will treat the problem as interactive image segmentation
      -------------------------------------------------------------------------------------

 
  







