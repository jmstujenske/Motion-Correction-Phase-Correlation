# Motion-Correction-Phase-Correlation
Standalone version of Matlab suite2p phase correlation motion correction, with subpixel registration


Includes a novel bidirectional correction algorithm which accounts for unevenness in bidirectional scanning offset along the scan axis. At least for Olympus 2p, this is quite helpful.


reg2P_standalone_fullstack_slowdrift is the core script that I would recommend for long recordings. It both handles motion correction from scanning artifacts and slow drift over minutes which can occur in some imaging sessions.
