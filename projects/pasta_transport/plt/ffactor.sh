gnuplot form_factors.gp
ffmpeg -r 5 -i form_factor%02d.png -pix_fmt yuv420p -vcodec libx264 out.avi
xdg-open out.avi
