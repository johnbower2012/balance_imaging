EXP_FOLDER='exp_run00/'
for((i=$1;i<$2;i++))
do
    fn=` printf "model_run%02d/" $i`
    echo Now running for directory ${fn}...
    cp -v calculations/data/${fn}I211_J211.dat calculations/balancefunction/${fn}prior/
    cp -v calculations/data/${fn}I2212_J2212.dat calculations/balancefunction/${fn}prior/
    cp -v calculations/data/${fn}I321_J2212.dat calculations/balancefunction/${fn}prior/
    cp -v calculations/data/${fn}I321_J321.dat calculations/balancefunction/${fn}prior/

    echo Now running for directory posterior/${fn}...
    echo "paste calculations/balancefunction/${fn}posterior/I211_J211.dat <(cut -d ' ' -f 2 calculations/data/exp_data/star_pipi_cut.dat) > calculations/balancefunction/${fn}posterior/pipi_exp.dat"
    paste calculations/balancefunction/${fn}posterior/I211_J211.dat <(cut -d ' ' -f 2 calculations/data/${EXP_FOLDER}star_pipi_cut.dat) > calculations/balancefunction/${fn}posterior/pipi_exp.dat
    echo "paste calculations/balancefunction/${fn}posterior/I2212_J2212.dat <(cut -d ' ' -f 2 calculations/data/exp_data/star_ppbar_cut.dat) > calculations/balancefunction/${fn}posterior/ppbar_exp.dat"
    paste calculations/balancefunction/${fn}posterior/I2212_J2212.dat <(cut -d ' ' -f 2 calculations/data/${EXP_FOLDER}star_ppbar_cut.dat) > calculations/balancefunction/${fn}posterior/ppbar_exp.dat
    echo "paste calculations/balancefunction/${fn}posterior/I321_J2212.dat <(cut -d ' ' -f 2 calculations/data/exp_data/star_pK_cut.dat) > calculations/balancefunction/${fn}posterior/pk_exp.dat"
    paste calculations/balancefunction/${fn}posterior/I321_J2212.dat <(cut -d ' ' -f 2 calculations/data/${EXP_FOLDER}star_pK_cut.dat) > calculations/balancefunction/${fn}posterior/pk_exp.dat
    echo "paste calculations/balancefunction/${fn}posterior/I321_J321.dat <(cut -d ' ' -f 2 calculations/data/exp_data/star_KK_cut.dat) > calculations/balancefunction/${fn}posterior/kk_exp.dat"
    paste calculations/balancefunction/${fn}posterior/I321_J321.dat <(cut -d ' ' -f 2 calculations/data/${EXP_FOLDER}star_KK_cut.dat) > calculations/balancefunction/${fn}posterior/kk_exp.dat


    echo Now running for directory prior/${fn}...
    echo "paste calculations/balancefunction/${fn}prior/I211_J211.dat <(cut -d ' ' -f 2 calculations/data/exp_data/star_pipi_cut.dat) > calculations/balancefunction/${fn}prior/pipi_exp.dat"
    paste calculations/balancefunction/${fn}prior/I211_J211.dat <(cut -d ' ' -f 2 calculations/data/${EXP_FOLDER}star_pipi_cut.dat) > calculations/balancefunction/${fn}prior/pipi_exp.dat
    echo "paste calculations/balancefunction/${fn}prior/I2212_J2212.dat <(cut -d ' ' -f 2 calculations/data/exp_data/star_ppbar_cut.dat) > calculations/balancefunction/${fn}prior/ppbar_exp.dat"
    paste calculations/balancefunction/${fn}prior/I2212_J2212.dat <(cut -d ' ' -f 2 calculations/data/${EXP_FOLDER}star_ppbar_cut.dat) > calculations/balancefunction/${fn}prior/ppbar_exp.dat
    echo "paste calculations/balancefunction/${fn}prior/I321_J2212.dat <(cut -d ' ' -f 2 calculations/data/exp_data/star_pK_cut.dat) > calculations/balancefunction/${fn}prior/pk_exp.dat"
    paste calculations/balancefunction/${fn}prior/I321_J2212.dat <(cut -d ' ' -f 2 calculations/data/${EXP_FOLDER}star_pK_cut.dat) > calculations/balancefunction/${fn}prior/pk_exp.dat
    echo "paste calculations/balancefunction/${fn}prior/I321_J321.dat <(cut -d ' ' -f 2 calculations/data/exp_data/star_KK_cut.dat) > calculations/balancefunction/${fn}prior/kk_exp.dat"
    paste calculations/balancefunction/${fn}prior/I321_J321.dat <(cut -d ' ' -f 2 calculations/data/${EXP_FOLDER}star_KK_cut.dat) > calculations/balancefunction/${fn}prior/kk_exp.dat

    echo "Ending directory ${fn}..."
    echo ""
done