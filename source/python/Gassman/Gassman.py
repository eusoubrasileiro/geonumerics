

def BulkModulus(Vp, Vs, constBulkDensity):
    """
    Calculate K modulus from Vp, Vs and bulk density of the rock
    Vp and Vs must be array's. (Remember units, be consistent)
    """
    
    return constBulkDensity*(Vp**2-4*Vs**2/3);


# to make fluid substituition, first we have to define:
# 1) porosity of rock
# 2) properties of the fluids (K_fl, Density_fl)
# 3) bulk modulus of mineral matrix K0
# 4) bulk modulus of porous rock frame K*



def OilRefDensityFromApi(Api):
    """
    Calculate the reference density from the API oil    classification
    (density units = g/cm^3)
    """
    return  141.5/(Api+131.5);


def OilDensityFromP(RefOilDensity, Pressure):
    """
    RefOilDensity is the oil density
    Calculate the oil density to an specific
    pressure (based on reference density api)
    units (pressure=MPa .. density=g/cm^3)
    McCain 1973 look
    """
    Density=RefOilDensity+(0.00277*Pressure-(Pressure**3)*1.71*(10**-7))*(RefOilDensity-1.15)**2;
    Density+=Pressure*3.49*(10**-4);

    return Density

#    kg/m^3 = 1000g/(100cm)^3= 1000g/10^5cm^3 = g/cm^3 * 10^-3
# pascal = Newton/m^2 = Kg*m/s^2 / m^2
# 1MPa equivalente = 1000*ton/m^2
# 1KPa 1ton/m^2

# 1g/cm^3 = (10^-3)kg/ (10^-2m)^3 =10^-3/10^-6= 10^3 kg/m^3 

def OilDensityFromT(DensityP, Temperature):
    """
    DensityP is the oil density
    Dodson and Standing (1945)    
    Density at a fixed pressure, expression varying with temperature
    (units density (g/cm^3) and Temperature (celcius))
    """
    Density=DensityP
    Density/=0.972+3.81*(10**-4)*(Temperature+17.78)**(1.175)
    return Density;


def OilDensityFromTP(RefDensity, Temperature, Pressure):
    """
    RefDensity is the oil density API
    Wrapper of OilDensityFromP and OilDensityFromT
    temperature in Celcius
    pressure in MPa
    """
    DensityP=OilDensityFromP(RefDensity, Pressure);
    return OilDensityFromT(DensityP, Temperature)

#Lithostatic pressure Density*gravity*H


def OilVpFromTP(OilDensity, Temperature, Pressure):
    """
    OilDensity is the oil density
    Oil velocity from density, temperature and pressure
    Vp in m/s
    temperature in Celcius
    pressure in MPa
    """
    T=Temperature;
    P=Pressure;
    rho=OilDensity;
    V = 2096*(rho/(2.6-rho))**0.5
    V += 4.64*P - 3.7*T
    V += 0.0115*( 4.12*(1.08*(rho**-1) -1)**0.5 - 1 )*T*P
    return V


def OilKFromTP(OilDensity, Temperature, Pressure):
    """
    Oil bulk modulus from temperature and pressure
    density in g/cm^3
    temperature in Celcius
    pressure in MPa
    uses OilVpFromTP 
    K in MPa
    """
    # P velocity only possible on fluids
    Vp = OilVpFromTP(OilDensity, Temperature, Pressure)
    # from bulk modulus definition.. (with Vs = 0)!!
    # converting density to Kg/m^3
    K = 1000*OilDensity*Vp**2
    # kg/m^3 * m^2 / s^2 = Pascal
    # convert to MPA divide by 1 million
    K /= 10**6     
    return K

def OilBulkModulusFromTP(RefOilDensity, Temperature, Pressure):
    """
    Oil bulk modulus from temperature and pressure
    density in g/cm^3 (API density)
    temperature in Celcius
    pressure in MPa
    uses OilVpFromTP 
    K in MPa
    """
    OilDensityAtTP=OilDensityFromTP(RefOilDensity, Temperature, Pressure)
    # P velocity only possible on fluids
    Vp = OilVpFromTP(OilDensityAtTP, Temperature, Pressure)
    # from bulk modulus definition.. (with Vs = 0)!!
    # converting density to Kg/m^3
    K = 1000*OilDensityAtTP*Vp**2
    # convert to MPA divide by 1 million
    K /= 10**6  
    return K

