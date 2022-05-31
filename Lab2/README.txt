Clara Llòria clara.lloria01@estudiant.upf.edu 218147
Ainhoa Pachón ainhoa.pachon01@estudiant.upf.edu 219371

Hemos intentado implementar el uso de GL_GREATER para el "volume lightning" en Deferred, pero aparecían unos artefactos que no hemos sabido solucionar y 
hemos dejado esa parte de código comentada.
Hemos encontrado un bug donde si quitamos las spotlights en deferred la ambient light se come el resto de luces. No hemos encontrado de donde sale el error.
También hemos visto que si te acercas mucho a algunas pointlights (no todas), la luz de esta desaparece.