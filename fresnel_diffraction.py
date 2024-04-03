import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.integrate import trapz

def menu_text():
    """This function prints a menu listing the various sections of the program which may be run"""
    print (30*"-", "MENU", 30*"-")
    print ("1. 1D Fresnel Diffraction")
    print ("2. 2D Fresnel Diffraction With Simpson's Rule")
    print ("3. 2D Fresnel Diffraction With Trapezoid Rule")
    print ("4. 1D Comparison: Simpson's and Trapezoid Rules")
    print ("5. Exit")
    print (67*"-")

def Intensity_function_1D(k,x,z,x1,x2,N):
    """This function calculates the intensity of light on a screen distance z from an 
    aperture width x1-x2 as a function of position on the screen x through integration
    using simpsons rule of N intervals"""
    Aperture_vals=np.linspace(x1, x2, N+1) #sets up array of aperture values with N intervals across the user defined limits
    h=(x2-x1)/N #calculates step size h
    Integral=0
    for i in range(len(Aperture_vals)): #cycles through aperture values
        if i==0 or i==N:
            a=(math.cos((k*((x-Aperture_vals[i])**2))/(2*z))+(1j*math.sin((k*((x-Aperture_vals[i])**2))/(2*z)))) #if first and last integrand entry, multiplied by 1
        elif i%2!=0:
            a=4*(math.cos((k*((x-Aperture_vals[i])**2))/(2*z))+(1j*math.sin((k*((x-Aperture_vals[i])**2))/(2*z)))) #if odd integrand entry, multiplied by 4
        else:
            a=2*(math.cos((k*((x-Aperture_vals[i])**2))/(2*z))+(1j*math.sin((k*((x-Aperture_vals[i])**2))/(2*z)))) #if even integrand entry that is not first and last, multiplied by 2
        Integral+=a
        E=(h/3)*Integral*(k/(2*math.pi*z)) #sum is multiplied by step size/3 to give the simpsons integration, then multiplied by a set of constants and variables to give the E field.
        Intensity=(8.85E-12)*(3E8)*((abs(E))**2) #modulus squared of the E field (same as complex times conjugate) multiplied by set of constants returns Intensity value
    return Intensity

def Intensity_function_x_2D(k,x,z,x1,x2,N):
    """This function performs an inbuilt simpsons rule integration across the aperture width x1'-x2'
    for each coordinate x of the screen"""
    Aperture_vals=np.linspace(x1, x2, N+1) #sets up array of aperture values with N intervals across the user defined limits
    Integrand=np.zeros(N+1,dtype=complex) #initiates complex array to accept values of the integrand (which are complex numbers)
    for i in range(N+1):
        Integrand[i]=np.exp(1j*(k/(2*z))*((x-Aperture_vals[i])**2)) #calculates value of the integrand due to each aperture coordinate
    Integral=simps(Integrand,Aperture_vals) #performs the simpsons integration with the aperture and integrand arrays
    return Integral

def Intensity_function_y_2D(k,y,z,y1,y2,N):
    """This function performs an inbuilt simpsons rule integration across the aperture width y1'-y2'
    for each coordinate y of the screen"""
    Aperture_vals=np.linspace(y1, y2, N+1) #sets up array of aperture values with N intervals across the user defined limits
    Integrand=np.zeros(N+1,dtype=complex) #initiates complex array to accept values of the integrand (which are complex numbers)
    for i in range(N+1):
        Integrand[i]=np.exp(1j*(k/(2*z))*((y-Aperture_vals[i])**2)) #calculates value of the integrand due to each aperture coordinate
    Integral=simps(Integrand,Aperture_vals) #performs the simpsons integration with the aperture and integrand arrays
    return Integral

def Intensity_function_x_2D_Trapezoid(k,x,z,x1,x2,N):
    """This function performs an inbuilt trapezoid rule integration across the aperture width x1'-x2'
    for each coordinate x of the screen"""
    Aperture_vals=np.linspace(x1, x2, N+1) #sets up array of aperture values with N intervals across the user defined limits
    Integrand=np.zeros(N+1,dtype=complex) #initiates complex array to accept values of the integrand (which are complex numbers)
    for i in range(N+1):
        Integrand[i]=np.exp(1j*(k/(2*z))*((x-Aperture_vals[i])**2)) #calculates value of the integrand due to each aperture coordinate
    Integral=trapz(Integrand,Aperture_vals) #performs the trapezoid integration with the aperture and integrand arrays
    return Integral

def Intensity_function_y_2D_Trapezoid(k,y,z,y1,y2,N):
    """This function performs an inbuilt trapezoid rule integration across the aperture width y1'-y2'
    for each coordinate y of the screen"""
    Aperture_vals=np.linspace(y1, y2, N+1) #sets up array of aperture values with N intervals across the user defined limits
    Integrand=np.zeros(N+1,dtype=complex) #initiates complex array to accept values of the integrand (which are complex numbers)
    for i in range(N+1):
        Integrand[i]=np.exp(1j*(k/(2*z))*((y-Aperture_vals[i])**2)) #calculates value of the integrand due to each aperture coordinate
    Integral=trapz(Integrand,Aperture_vals) #performs the trapezoid integration with the aperture and integrand arrays
    return Integral

def input_checker(a):                                                 
    """This function checks whether an input is a float"""
    try: 
        float(a) #attempt to take the float of argument a, if successful, returns a True
        return True
    except ValueError: #any input leading to a ValueError, e.g. a string, returns a False
        return False
    
def input_positive_checker(a): #used for the wavelength and distance variables which should be positive values                                   
    """This function checks whether an input is a positive float"""
    try:
        float(a)
        if float(a)<0: #additional condition to the above function that any negative inputs returns a False
            return False
        return True
    except ValueError:
        return False    
    
def N_check(a): #used for the intervals N which must be a positive integer
    """This function checks whether an input is a positive integer>0."""
    try:
        int(a) #repeat of above function with integers tried instead of floats                                                           
        if int(a)<0:                                                          
            return False                                                       
        return True                                                            
    except ValueError:
        return False
    
Loop=True #loop for menu selection initiated
print("Welcome to the program!")
while Loop:
    print ()
    menu_text() #menu detailing the various menu selections. Defined as a function above
    Menu_selection=input("Please enter one of the above Menu Options [1-5]: ")
    print()
    
    if Menu_selection=="1":
        wavelength=input("Welcome\nPlease insert a value for wavelength, λ, in metres:")
        while input_positive_checker(wavelength)!=True:
            wavelength=input("I'm sorry, this is not a valid wavelength\nPlease insert a positive wavelength in metres:")
        k=(2*math.pi)/float(wavelength) #checks user input to ensure positive float inputted, only after this is satisfied will the float of the input be taken
        x1=input("Thank you\nPlease insert the starting coordinate of the aperture, x1', in metres:")
        while input_checker(x1)!=True:
            x1=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid starting coordinate of the aperture in metres:")
        x1=float(x1)
        x2=input("Thank you\nPlease insert the ending coordinate of the aperture, x2', in metres:")
        while input_checker(x2)!=True:
            x2=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid ending coordinate of the aperture in metres:")
        x2=float(x2)
        z=input("Thank you\nPlease insert the screen distance to the aperture, z, in metres:")
        while input_positive_checker(z)!=True:
            z=input("I'm sorry, this is not a valid distance\nPlease insert a valid positive screen distance to the aperture in metres:")
        z=float(z)
        N=input("Please insert an even interval integer for integration, N:")
        while N_check(N)!=True:
            N=input("I'm sorry, this is not a valid interval\nPlease insert a valid even integer interval for integration:")
        N=int(N)
        xmin=input("Thank you\nPlease insert the desired minimum screen coordinate for the intensity plot, xmin, in metres:")
        while input_checker(xmin)!=True:
            xmin=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid minimum screen coordinate in metres:")
        xmin=float(xmin)
        xmax=input("Thank you\nPlease insert the desired maximum screen coordinate for the intensity plot, xmax, in metres:")
        while input_checker(xmax)!=True:
            xmax=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid maximum screen coordinate in metres:")
        xmax=float(xmax)
        NumPoints=200
        xvals=np.linspace(xmin,xmax,NumPoints) #sets up an array of 200 values from xmin to xmax of the user defined screen limits
        yvals=np.zeros(NumPoints)
        for i in range(NumPoints):
            yvals[i]=Intensity_function_1D(k,xvals[i],z,x1,x2,N) #calculates intensity for each screen coordinate in the array using self built simpsons integration function
        plt.plot(xvals,yvals)
        plt.xlabel("Screen coordinate (m)")
        plt.ylabel("Relative Intensity")
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0)) #sets both axes ticks to scientific notation for improved formatting
        plt.show()
                
    elif Menu_selection=="2":
        wavelength=input("Welcome\nPlease insert a value for wavelength, λ, in metres:")
        while input_positive_checker(wavelength)!=True:
            wavelength=input("I'm sorry, this is not a valid wavelength\nPlease insert a positive wavelength in metres:")
        k=(2*math.pi)/float(wavelength) #As before, checks user input to ensure positive float inputted, only after this is satisfied will the float of the input be taken
        x1=input("Thank you\nPlease insert the starting x coordinate of the aperture, x1', in metres:")
        while input_checker(x1)!=True:
            x1=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid starting x coordinate of the aperture in metres:")
        x1=float(x1)
        x2=input("Thank you\nPlease insert the ending x coordinate of the aperture, x2', in metres:")
        while input_checker(x2)!=True:
            x2=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid ending x coordinate of the aperture in metres:")
        x2=float(x2)
        y1=input("Thank you\nPlease insert the starting y coordinate of the aperture, y1', in metres:")
        while input_checker(y1)!=True:
            y1=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid starting y coordinate of the aperture in metres:")
        y1=float(y1)
        y2=input("Thank you\nPlease insert the ending y coordinate of the aperture, y2', in metres:")
        while input_checker(y2)!=True:
            y2=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid ending y coordinate of the aperture in metres:")
        y2=float(y2)
        z=input("Thank you\nPlease insert the screen distance to the aperture, z, in metres:")
        while input_positive_checker(z)!=True:
            z=input("I'm sorry, this is not a valid distance\nPlease insert a valid positive screen distance to the aperture in metres:")
        z=float(z)
        N=input("Please insert an even interval integer for integration, N:")
        while N_check(N)!=True:
            N=input("I'm sorry, this is not a valid interval\nPlease insert a valid even integer interval for integration:")
        N=int(N)
        NumPoints=200
        xmin=input("Thank you\nPlease insert the desired minimum x screen coordinate for the intensity plot, xmin, in metres:")
        while input_checker(xmin)!=True:
            xmin=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid minimum x screen coordinate in metres:")
        xmin=float(xmin)
        xmax=input("Thank you\nPlease insert the desired maximum x screen coordinate for the intensity plot, xmax, in metres:")
        while input_checker(xmax)!=True:
            xmax=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid maximum x screen coordinate in metres:")
        xmax=float(xmax)
        ymin=input("Thank you\nPlease insert the desired minimum y screen coordinate for the intensity plot, ymin, in metres:")
        while input_checker(ymin)!=True:
            ymin=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid minimum y screen coordinate in metres:")
        ymin=float(ymin)
        ymax=input("Thank you\nPlease insert the desired maximum y screen coordinate for the intensity plot, ymax, in metres:")
        while input_checker(ymax)!=True:
            ymax=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid maximum y screen coordinate in metres:")
        ymax=float(ymax)
        xvals=np.linspace(xmin,xmax,100)
        yvals=np.linspace(ymin,ymax,100) #x screen coordiante, y screen coordinate arrays generated from the user defined screen limits
        NumPoints=100
        zvals=np.zeros((NumPoints,NumPoints)) 
        for i in range(len(xvals)):
            for j in range(len(yvals)): #for each x coordinate, all the intensity values at every y coordinate calculated. This is then continued cycling through each x coordinate.
                zvals[i,j]=8.85E-12*3E8*(abs((k/(2*np.pi*z))*Intensity_function_x_2D(k,xvals[i],z,x1,x2,N)*Intensity_function_y_2D(k,yvals[j],z,y1,y2,N)))**2
        plt.imshow(zvals, extent=[xmin,xmax,ymin,ymax]) #plots image to the extent of the user defined screen limits
        plt.locator_params(axis="both", nbins=5) #sets both axes to contain only 5 ticks
        cbar=plt.colorbar()
        cbar.set_label("Relative Intensity", rotation=270, labelpad=15) #labels colourbar. Label is rotated to vertical and set a suitable distance to the colourbar itself
        cbar.formatter.set_powerlimits((0, 0)) #Colourbar given scientific notation for improved formatting 
        cbar.update_ticks()
        plt.xlabel("x screen coordinate (m)")
        plt.ylabel("y screen coordinate (m)")
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0)) #both axes given scientific notation for improved formatting
        plt.show()
            
    elif Menu_selection=="3":
        wavelength=input("Welcome\nPlease insert a value for wavelength, λ, in metres:")
        while input_positive_checker(wavelength)!=True:
            wavelength=input("I'm sorry, this is not a valid wavelength\nPlease insert a positive wavelength in metres:")
        k=(2*math.pi)/float(wavelength) #As before, checks user input to ensure positive float inputted, only after this is satisfied will the float of the input be taken.
        x1=input("Thank you\nPlease insert the starting x coordinate of the aperture, x1', in metres:")
        while input_checker(x1)!=True:
            x1=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid starting x coordinate of the aperture in metres:")
        x1=float(x1)
        x2=input("Thank you\nPlease insert the ending x coordinate of the aperture, x2', in metres:")
        while input_checker(x2)!=True:
            x2=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid ending x coordinate of the aperture in metres:")
        x2=float(x2)
        y1=input("Thank you\nPlease insert the starting y coordinate of the aperture, y1', in metres:")
        while input_checker(y1)!=True:
            y1=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid starting y coordinate of the aperture in metres:")
        y1=float(y1)
        y2=input("Thank you\nPlease insert the ending y coordinate of the aperture, y2', in metres:")
        while input_checker(y2)!=True:
            y2=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid ending y coordinate of the aperture in metres:")
        y2=float(y2)
        z=input("Thank you\nPlease insert the screen distance to the aperture, z, in metres:")
        while input_positive_checker(z)!=True:
            z=input("I'm sorry, this is not a valid distance\nPlease insert a valid positive screen distance to the aperture in metres:")
        z=float(z)
        N=input("Please insert an even interval integer for integration, N:")
        while N_check(N)!=True:
            N=input("I'm sorry, this is not a valid interval\nPlease insert a valid even integer interval for integration:")
        N=int(N)
        NumPoints=200
        xmin=input("Thank you\nPlease insert the desired minimum x screen coordinate for the intensity plot, xmin, in metres:")
        while input_checker(xmin)!=True:
            xmin=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid minimum x screen coordinate in metres:")
        xmin=float(xmin)
        xmax=input("Thank you\nPlease insert the desired maximum x screen coordinate for the intensity plot, xmax, in metres:")
        while input_checker(xmax)!=True:
            xmax=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid maximum x screen coordinate in metres:")
        xmax=float(xmax)
        ymin=input("Thank you\nPlease insert the desired minimum y screen coordinate for the intensity plot, ymin, in metres:")
        while input_checker(ymin)!=True:
            ymin=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid minimum y screen coordinate in metres:")
        ymin=float(ymin)
        ymax=input("Thank you\nPlease insert the desired maximum y screen coordinate for the intensity plot, ymax, in metres:")
        while input_checker(ymax)!=True:
            ymax=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid maximum y screen coordinate in metres:")
        ymax=float(ymax)
        xvals=np.linspace(xmin,xmax,100)
        yvals=np.linspace(ymin,ymax,100) #x screen coordiante, y screen coordinate arrays generated from the user defined screen limits
        NumPoints=100
        zvals=np.zeros((NumPoints,NumPoints))
        for i in range(len(xvals)):
            for j in range(len(yvals)): #for each x coordinate, all the intensity values at every y coordinate calculated. This is then continued cycling through each x coordinate.
                zvals[i,j]=8.85E-12*3E8*(abs((k/(2*np.pi*z))*Intensity_function_x_2D_Trapezoid(k,xvals[i],z,x1,x2,N)*Intensity_function_y_2D_Trapezoid(k,yvals[j],z,y1,y2,N)))**2
        plt.imshow(zvals, extent=[xmin,xmax,ymin,ymax]) #plots image to the extent of the user defined screen limits
        plt.locator_params(axis="both", nbins=5) #sets both axes to contain only 5 ticks
        cbar=plt.colorbar()
        cbar.set_label("Relative Intensity", rotation=270, labelpad=15) #labels colourbar. Label is rotated to vertical and set a suitable distance to the colourbar itself
        cbar.formatter.set_powerlimits((0, 0)) #Colourbar given scientific notation for improved formatting
        cbar.update_ticks()
        plt.xlabel("x screen coordinate (m)")
        plt.ylabel("y screen coordinate (m)")
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0)) #both axes given scientific notation for improved formatting
        plt.show()
        
    elif Menu_selection=="4":
        wavelength=input("Welcome\nPlease insert a value for wavelength, λ, in metres:")
        while input_positive_checker(wavelength)!=True:
            wavelength=input("I'm sorry, this is not a valid wavelength\nPlease insert a positive wavelength in metres:")
        k=(2*math.pi)/float(wavelength) #As before, checks user input to ensure positive float inputted, only after this is satisfied will the float of the input be taken.
        x1=input("Thank you\nPlease insert the starting coordinate of the aperture, x1', in metres:")
        while input_checker(x1)!=True:
            x1=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid starting coordinate of the aperture in metres:")
        x1=float(x1)
        x2=input("Thank you\nPlease insert the ending coordinate of the aperture, x2', in metres:")
        while input_checker(x2)!=True:
            x2=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid ending coordinate of the aperture in metres:")
        x2=float(x2)
        z=input("Thank you\nPlease insert the screen distance to the aperture, z, in metres:")
        while input_positive_checker(z)!=True:
            z=input("I'm sorry, this is not a valid distance\nPlease insert a valid positive screen distance to the aperture in metres:")
        z=float(z)
        N=input("Please insert an even interval integer for integration, N:")
        while N_check(N)!=True:
            N=input("I'm sorry, this is not a valid interval\nPlease insert a valid even integer interval for integration:")
        N=int(N)
        NumPoints=200
        xmin=input("Thank you\nPlease insert the desired minimum screen coordinate for the intensity plot, xmin, in metres:")
        while input_checker(xmin)!=True:
            xmin=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid minimum screen coordinate in metres:")
        xmin=float(xmin)
        xmax=input("Thank you\nPlease insert the desired maximum screen coordinate for the intensity plot, xmax, in metres:")
        while input_checker(xmax)!=True:
            xmax=input("I'm sorry, this is not a valid coordinate\nPlease insert a valid maximum screen coordinate in metres:")
        xmax=float(xmax)
        xvals=np.linspace(xmin,xmax,NumPoints) #sets up an array of 200 values from xmin to xmax of the user defined screen limits
        yvals_simps=np.zeros(NumPoints)
        yvals_trapz=np.zeros(NumPoints)
        ydiff=np.zeros(NumPoints) #initiates arrays for simpsons, trapezoid and the difference
        for i in range(NumPoints):
            yvals_simps[i]=8.85E-12*3E8*(abs((k/(2*np.pi*z))*Intensity_function_x_2D(k,xvals[i],z,x1,x2,N)))**2 #calculates intensity for each screen coordinate in the array using inbuilt simpsons integration
            yvals_trapz[i]=8.85E-12*3E8*(abs((k/(2*np.pi*z))*Intensity_function_x_2D_Trapezoid(k,xvals[i],z,x1,x2,N)))**2 #calculates intensity for each screen coordinate in the array using inbuilt trapezoid integration
            ydiff[i]=abs(yvals_simps[i]-yvals_trapz[i]) #calculates difference between values produced from simpson's and trapezoid integrations
        plt.plot(xvals,yvals_simps, label="Simpsons")
        plt.plot(xvals,yvals_trapz, label="Trapezoid")
        plt.xlabel("Screen coordinate (m)")
        plt.ylabel("Relative Intensity")
        plt.legend()
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0)) #both axes given scientific notation for improved formatting
        plt.show()
        
        plt.plot(xvals, ydiff)
        plt.xlabel("Screen coordinate (m)")
        plt.ylabel("Relative Intensity")
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0)) #both axes given scientific notation for improved formatting
        plt.show()
        
    elif Menu_selection=="5":
        print("Thank you for using this program")
        break #ends the menu selection loop which in turn ends the program
    
    else: #accepts arguments which are not [1-5] so that the menu selection loop is not accidentally broken
        print("Apologies but I do not recognise this selection")

