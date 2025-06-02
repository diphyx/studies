import torch

def score_binding_site():
    if not torch.cuda.is_available():
        print("CUDA is not available. Please check your environment.")
        return

    device = torch.device("cuda")

    # Simulate a protein pocket as a 3D grid with "charge" intensities
    protein_grid = torch.zeros((32, 32, 32), device=device)
    protein_grid[15:18, 15:18, 15:18] = torch.tensor(1.0, device=device)  # binding pocket

    # Simulate a ligand field (probe) with complementary shape
    ligand_probe = torch.ones((3, 3, 3), device=device)

    # Perform cross-correlation (dot product) to find best binding location
    score_map = torch.nn.functional.conv3d(
        protein_grid.unsqueeze(0).unsqueeze(0),  # [N, C, D, H, W]
        ligand_probe.unsqueeze(0).unsqueeze(0),  # [out_channels, in_channels, kD, kH, kW]
        padding=1
    ).squeeze()

    best_score = score_map.max().item()
    best_position = torch.nonzero(score_map == score_map.max(), as_tuple=False)

    print(f"Best binding score: {best_score:.2f}")
    print(f"Best binding site grid position: {best_position.tolist()}")

if __name__ == "__main__":
    score_binding_site()
