

class KMC:
    def __init__(self, system, TDlib):
        self.system = system
        self.TDlib = TDlib
        self.missingTDs = np.array([]) # Add the missing fingerprints to this list
    
    def run(n):
        """
        Simulates n electrons, and return True or False depending on whether any TD values were missing. If no TD value was missing, return False, otherwise return True.
        """
        tF = 0
        for i in range(n):
            tF += self.simulate_electron() # simulate_electron() returns False (0) of no TD value was missing, and True (1) if a value was missing
        
        return bool(tF)
        
    
    def simulate_electron():
        ### TEMPORARILY A VIRTUAL FUNCTION ###
        raise NotImplementedError()
        
        #
        # 1. Determine where to hit structure
        # 2. Do something with the cross-section of atoms nearby (perhaps a search?), to determine which atom(s) are affected
        # 2.1 Check whether a TD value exists for this/these atom(s). If not, return 1
        # 3. Determine how much energy is transferred, if energy is transferred
        # - Use mr. PHD's distribution for this
        # 4. If above threshold value, remove the atom
        # 4.1 Maybe implement some other behaviour, such as removing neighboring atoms
        # 4.2 Also perhaps add the energy to the system in case it does not exceed the TD
        #
        # Return 0 at the end (if the TD value for the atom exists)
        #
    
    def cross_section():
        ### TEMPORARILY A VIRTUAL FUNCTION ###
        raise NotImplementedError()
        
    def cross_dummy():
        """
        Dummy-function to imitate calculating the interaction cross-section of an atom.
        """
    
    def determine_hit():
        """
        Determines where an electron enters the structures, and determines whether it'll intersect any cross-section for an atom. It then returns the ID of this atom if True, otherwise it returns False.
        """
        
    def get_TD():
        """
        Using an id, calculate the fingerprint of the atom and return its TD value, using the TD library. If no TD value exists for the given fingerprint, log this and return False.
        """
        ### TEMPORARILY A VIRTUAL FUNCTION ###
        raise NotImplementedError()
    
    def get_missing_TDs():
        """
        Returns all the fingerprints with missing TD values.
        """
        return self.missingTDs
    
    
