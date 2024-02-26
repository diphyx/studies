import org.openscience.cdk.Atom;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.config.Elements;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import java.io.FileWriter; // Import the FileWriter class

public class CDKTest {
    public static void main(String[] args) throws Exception {
        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        IAtomContainer molecule = builder.newInstance(IAtomContainer.class);
        
        IAtom carbon = builder.newInstance(IAtom.class, Elements.CARBON);
        molecule.addAtom(carbon);

        // Ensure the molecule has all atom types perceived
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);

        // Add implicit hydrogens.
        CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
        adder.addImplicitHydrogens(molecule);
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);

        // Generating SMILES
        SmilesGenerator sg = SmilesGenerator.generic();
        String smiles = sg.create(molecule);

        // Write the SMILES to a file
        try (FileWriter writer = new FileWriter("smiles_output.txt")) {
            writer.write("SMILES: " + smiles);
            System.out.println("Successfully wrote SMILES to file.");
        } catch (Exception e) {
            System.out.println("An error occurred while writing to the file.");
            e.printStackTrace();
        }
    }
}
