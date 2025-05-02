Qav_mse = mean((Qav - QavQP)^2)
Qmt_mse = mean((Qmt - QmtQP)^2)
Qpv_mse = mean((Qpv - QpvQP)^2)
Qtc_mse = mean((Qtc - QtcQP)^2)
Vao_mse = mean((Vao - VaoQP)^2)
Vlv_mse = mean((Vlv - VlvQP)^2)
Vpa_mse = mean((Vpa - VpaQP)^2)
Vpu_mse = mean((Vpu - VpuQP)^2)
Vrv_mse = mean((Vrv - VrvQP)^2)
Vvc_mse = mean((Vvc - VvcQP)^2)
Vspt_mse = mean((Vspt - VsptQP)^2)

Meanmse = data.frame(Qav_mse, Qmt_mse, Qpv_mse,
                 Qtc_mse, Vao_mse, Vlv_mse,
                 Vpa_mse, Vpu_mse, Vrv_mse,
                 Vvc_mse, Vspt_mse)
View(Meanmse)

Qav_msep3 = mean((Qav - QavQPp3)^2)
Qmt_msep3 = mean((Qmt - QmtQPp3)^2)
Qpv_msep3 = mean((Qpv - QpvQPp3)^2)
Qtc_msep3 = mean((Qtc - QtcQPp3)^2)
Vao_msep3 = mean((Vao - VaoQPp3)^2)
Vlv_msep3 = mean((Vlv - VlvQPp3)^2)
Vpa_msep3 = mean((Vpa - VpaQPp3)^2)
Vpu_msep3 = mean((Vpu - VpuQPp3)^2)
Vrv_msep3 = mean((Vrv - VrvQPp3)^2)
Vvc_msep3 = mean((Vvc - VvcQPp3)^2)
Vspt_msep3 = mean((Vspt - VsptQPp3)^2)

p3mse = data.frame(Qav_msep3, Qmt_msep3, Qpv_msep3,
                     Qtc_msep3, Vao_msep3, Vlv_msep3,
                     Vpa_msep3, Vpu_msep3, Vrv_msep3,
                     Vvc_msep3, Vspt_msep3)

View(p3mse)

Qav_msem3 = mean((Qav - QavQPm3)^2)
Qmt_msem3 = mean((Qmt - QmtQPm3)^2)
Qpv_msem3 = mean((Qpv - QpvQPm3)^2)
Qtc_msem3 = mean((Qtc - QtcQPm3)^2)
Vao_msem3 = mean((Vao - VaoQPm3)^2)
Vlv_msem3 = mean((Vlv - VlvQPm3)^2)
Vpa_msem3 = mean((Vpa - VpaQPm3)^2)
Vpu_msem3 = mean((Vpu - VpuQPm3)^2)
Vrv_msem3 = mean((Vrv - VrvQPm3)^2)
Vvc_msem3 = mean((Vvc - VvcQPm3)^2)
Vspt_msem3 = mean((Vspt - VsptQPm3)^2)

m3mse = data.frame(Qav_msem3, Qmt_msem3, Qpv_msem3,
                   Qtc_msem3, Vao_msem3, Vlv_msem3,
                   Vpa_msem3, Vpu_msem3, Vrv_msem3,
                   Vvc_msem3, Vspt_msem3)

View(m3mse)
