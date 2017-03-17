source ~/spinnaker_setup.sh
make clean
make || exit $?

cd ../IHC_AN_MSR_float

source ~/spinnaker_setup.sh
make clean
make || exit $?

cd ../IHC_AN_HSR_float

source ~/spinnaker_setup.sh
make clean
make || exit $?
